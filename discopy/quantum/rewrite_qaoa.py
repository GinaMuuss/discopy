from pyzx import VertexType
from fractions import Fraction
from pyzx.graph.base import VT
from sympy import sympify
from typing import List, Tuple
import pyzx as zx
from discopy.monoidal import Diagram, LocalSum
import sympy
import discopy
from discopy.quantum.zx import Diagram as discoZxDiag
from discopy.quantum.zx import X, Z, Id, Scalar, SWAP, PRO


def match_center_pi(graph, symbol) -> List[Tuple[VT, VT]]:
    """
    This function searches for Z spiders with a pi phase which have a symbol X spider to the right we can commute through.
    This is not the best subgraph search (actually it is rather bad....)
    Also: this could output also the part on the right

    Parameters
    ----------
    graph:
        pyzx graph
    symbol:
        sympy symbol which we match for

    Returns
    -------
    a List of candidates as tuples (candidate, neighbour)
    """
    phases = graph.phases()
    types = graph.types()
    c = []
    for v in graph.vertices():
        if phases[v] == Fraction(1, 1) and types[v] == zx.VertexType.Z:
            for n in graph.neighbors(v):
                if symbol in sympify(phases[n]).free_symbols:
                    c.append((v, n))
                    break
    return c


def permutate_pi_through(graph, candidates: List[Tuple[VT, VT]]):
    """
    this will permutate pi spiders through the symbol spiders

    Parameters
    ----------
    graph:
        pyzx graph
    candidates: List[Tuple[VT, VT]]
        a List of candidates as tuples (candidate, neighbour)

    Returns
    -------
    graph
    """
    phases = graph.phases()
    types = graph.types()
    print("candidates", candidates)
    edges_to_add = []
    edges_to_remove = []
    for c, n in candidates:
        # prepre0 -> pre0 -> r -> pre1
        # prepre0 -> r -> pre0 -> pre1
        other_neigh_c = list(graph.neighbors(c))
        other_neigh_c.remove(n)
        other_neigh_n = list(graph.neighbors(n))
        other_neigh_n.remove(c)

        for other_nei in other_neigh_c:
            edges_to_remove.append((c, other_nei))
            edges_to_add.append((n, other_nei))
        for other_nei in other_neigh_n:
            edges_to_remove.append((n, other_nei))
            edges_to_add.append((c, other_nei))
        graph.set_phase(n, -1 * phases[n])

    print("edges_to_add", edges_to_add)
    print("edges_to_remove", edges_to_remove)
    for e in edges_to_remove:
        graph.remove_edge(e)
    for e in edges_to_add:
        graph.add_edge(e)
    return graph


def find_candidates_for_pull_symbol_to_scalar(
    graph: discopy.quantum.zx.Diagram, inner_symbol
):
    """
    at the moment just returns all X spiders with the given symbol
    """
    candidates = []
    for i, box in enumerate(graph.boxes):
        if (
            isinstance(box, X)
            and hasattr(box, "phase")
            and hasattr(sympify(box.phase), "free_symbols")
            and inner_symbol in sympify(box.phase).free_symbols
        ):
            candidates.append(i)
    return candidates


def sum_from_Xspider(box):
    terms_left = Id(1) @ Scalar(sympy.cos(box.phase / 2))
    terms_right = X(1, 1, 0.5) @ Scalar(sympy.I * sympy.sin(box.phase / 2))
    return LocalSum([terms_left, terms_right])


def combine_sums(diag, sum_index_first):
    """This does not make even the most basic sanity checks!
    Parameters
    ----------
    diag: Diagram
        sum_index_first: int
        the index of the first sum, the second sum has to be directly behind
        the one behind is ASSUMED to be kroneker and not composed
    """
    assert isinstance(diag.boxes[sum_index_first], LocalSum)
    assert isinstance(diag.boxes[sum_index_first + 1], LocalSum)
    terms = []
    for term_left in diag.boxes[sum_index_first].terms:
        for term_right in diag.boxes[sum_index_first + 1].terms:
            terms += [term_left @ term_right]
    l = LocalSum(terms)
    new_boxes = diag.boxes[:sum_index_first] + [l] + diag.boxes[sum_index_first + 2 :]
    return discoZxDiag(
        diag.dom,
        diag.cod,
        new_boxes,
        diag.offsets[: sum_index_first + 1] + diag.offsets[sum_index_first + 2 :],
    )


def find_combinable_sums(diag):
    candidates = []
    for i, box in enumerate(diag.boxes[:-1]):
        if isinstance(box, LocalSum) and isinstance(diag.boxes[i + 1], LocalSum):
            candidates.append(i)
    return candidates


def simplify_inner(diagram: discopy.quantum.zx.Diagram, inner_symbol, outer_symbol):
    #diagram.draw()
    pyzx_final = diagram.to_pyzx()
    zx.simplify.clifford_simp(pyzx_final)
    #for i in range(5):
        #zx.simplify.spider_simp(pyzx_final)
        #zx.simplify.id_simp(pyzx_final)
    d = discopy.quantum.zx.Diagram.from_pyzx(pyzx_final)
    d.draw()
    return d


def simplify_qaoa(diagram: discopy.quantum.zx.Diagram, inner_symbol, outer_symbol):
    pyzx_final = diagram.to_pyzx()

    # Look for the center pi spiders
    # TODO: will this always be pi spiders? how can the center look?
    candidates = match_center_pi(pyzx_final, inner_symbol)
    pyzx_final = permutate_pi_through(pyzx_final, candidates)
    pyzx_final.normalize()

    # Do some basic simpification
    zx.simplify.spider_simp(pyzx_final)
    zx.simplify.id_simp(pyzx_final)
    disco_diag = discopy.quantum.zx.Diagram.from_pyzx(pyzx_final)

    # find some spiders with symbols where we pull the symbol into a scalar
    candidates = find_candidates_for_pull_symbol_to_scalar(disco_diag, inner_symbol)
    new_boxes = disco_diag.boxes
    for candidate in candidates:
        new_boxes[candidate] = sum_from_Xspider(disco_diag.boxes[candidate])
    new_diag = discoZxDiag(
        disco_diag.dom, disco_diag.cod, new_boxes, disco_diag.offsets
    )

    # we combine adjacent sums
    candidates = find_combinable_sums(new_diag)
    while len(candidates) > 0:
        new_diag = combine_sums(new_diag, candidates[0])
        candidates = find_combinable_sums(new_diag)

    # distribute the sum outward and simplify independetly
    candidates = []
    for i, box in enumerate(new_diag.boxes[:-1]):
        if isinstance(box, LocalSum):
            candidates.append(i)
    # TODO: handle more than one sum
    assert len(candidates) == 1
    candidate = candidates[0]
    results = []
    for term in new_diag.boxes[candidate].terms:
        dist_diag = discoZxDiag(
            new_diag.dom,
            new_diag.cod,
            new_diag.boxes[:candidate] + [term] + new_diag.boxes[candidate + 1 :],
            new_diag.offsets,
        )
        dist_diag = dist_diag.upgrade(
            discopy.quantum.zx.Functor(
                lambda x: x, lambda f: f, ob_factory=PRO, ar_factory=discoZxDiag
            )(dist_diag)
        )
        d = simplify_inner(dist_diag, inner_symbol, outer_symbol)
        results.append(d)

    return LocalSum(results)
