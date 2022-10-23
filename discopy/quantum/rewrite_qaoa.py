import itertools
from pyzx import VertexType, draw
import pyzx
from pyzx.graph.graph_s import GraphS
from fractions import Fraction
from pyzx.graph.base import VT
from pyzx.utils import EdgeType
from sympy import sympify
from typing import List, Tuple
import pyzx as zx
from discopy.monoidal import Diagram, LocalSum
import sympy
import discopy
from discopy.quantum.zx import Diagram as discoZxDiag
from discopy.quantum.zx import X, Z, Id, Scalar, SWAP, PRO
from queue import Queue


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


def find_cycles(diagram):
    # this is the worst cicle finder ever...
    # TODO: use a better algorithm to find 4-cycles
    # this only finds one cycle...
    unvisited_nodes = set(diagram.vertices())
    previous = {x: None for x in diagram.vertices()}
    stack = []  # tuples of nextnode and where we are comming from
    while len(stack) > 0 or len(unvisited_nodes) > 0:
        if len(stack) == 0:
            a = unvisited_nodes.pop()
            stack.append((a, None))
            unvisited_nodes.add(a)

        current, vfrom = stack.pop()
        for nei in diagram.neighbors(current):
            if nei == vfrom:
                continue
            stack.append((nei, current))

        if vfrom:
            previous[current] = previous[vfrom] + [vfrom]
        else:
            previous[current] = []

        if current in unvisited_nodes:
            unvisited_nodes.remove(current)
        else:
            return previous[current]


def find_bialg_reverse(diagram: GraphS):
    """
    Atm this only finds one 4 cycle, maybe not the best thing, but should work for now
    """
    cycle = list(find_cycles(diagram))
    for i, v in enumerate(cycle):
        a = diagram.edge_type(
            diagram.edge(v, cycle[i + 1 if i < len(cycle) - 1 else 0])
        )
        assert a == EdgeType.HADAMARD
    # TODO: sanity check that this 4 cycle is actually relevant...
    return cycle


def replace_bialg_reverse(diagram: GraphS, cycle: List[VertexType]):
    # TODO: double check that X spiders
    assert len(cycle) == 4

    # the first step ist to unspider, s.t. we can apply bialg
    for v in cycle:
        if diagram.phase(v) != 0:
            pyzx.rules.unspider(diagram, [v, []])
    
    print("unspidered the thing")
    draw(diagram)

    # first we collect the connections for the even and odd ones
    connect_to_even = []
    connect_to_odd = []
    # I know connect_to_even and connect_to_odd seem fliped, but the bialgebra rule does flip the assignment!
    for i in range(0, len(cycle), 2):
        connect_to_odd += [
            (x, diagram.edge_type(diagram.edge(x, cycle[i])))
            for x in diagram.neighbors(cycle[i])
            if x not in cycle
        ]
    for i in range(1, len(cycle), 2):
        connect_to_even += [
            (x, diagram.edge_type(diagram.edge(x, cycle[i])))
            for x in diagram.neighbors(cycle[i])
            if x not in cycle
        ]

    # now we remove the vertices
    diagram.remove_vertices(cycle)

    # now we add two now ones and their respective edges
    vEven, vOdd = diagram.add_vertices(2)
    diagram.set_type(vEven, VertexType.X)
    diagram.set_type(vOdd, VertexType.Z)
    diagram.add_edges([diagram.edge(vOdd, vEven)], EdgeType.SIMPLE)

    for v2, eType in connect_to_odd:
        realtype = EdgeType.HADAMARD if eType == EdgeType.SIMPLE else EdgeType.SIMPLE
        diagram.add_edges([diagram.edge(vOdd, v2)], realtype)
    for v2, eType in connect_to_even:
        diagram.add_edges([diagram.edge(vEven, v2)], eType)
    return diagram


def bad_heuristic_are_phases_zero(diagram: GraphS, symbol, step_size=0.1):
    """
    This is a bad heuristic.
    It will go through all phases and enter values for the inner_symbol from 0 to 2 with step_size.
    0 to 2 in pyzx means 0 to 2*pi.
    If the phase is 0 for all values entered it will reset it to said value.
    This is problamatic for several reasons and not even that fast.
    """
    phases = diagram.phases()
    for v in diagram.vertices():
        phase = sympy.Mod(phases[v], 2)
        print("vertex", v, phases[v], phase, phase.free_symbols)
        if symbol not in phase.free_symbols:
            #print("reset phase for vertex, because not in free_symbols", v, "to", phase, type(phase),  isinstance(phase, sympy.core.numbers.Number))
            if isinstance(phase, sympy.core.numbers.Number):
                diagram.set_phase(v, Fraction(float(phase)))
            continue
        all_zero = True
        val = 0
        base_val = phase.subs({symbol: val}) 
        while val <= 2:
            print(phase, val, phase.subs({symbol: val}))
            if phase.subs({symbol: val}) != base_val:
                all_zero = False
                break
            val += abs(step_size)
        if all_zero:
            #print("reset phase for vertex", v, "to", base_val)
            diagram.set_phase(v, base_val)
    return diagram


def simplify_inner(diagram: discopy.quantum.zx.Diagram, inner_symbol, outer_symbol):
    # diagram.draw()
    pyzx_final = diagram.to_pyzx()
    print("Before doing anything")
    draw(pyzx_final)
    zx.simplify.clifford_simp(pyzx_final)
    print("Before the bialg replace")
    draw(pyzx_final)
    cycle = find_bialg_reverse(pyzx_final)
    pyzx_final = replace_bialg_reverse(pyzx_final, cycle)
    print("After the bialg replace")
    draw(pyzx_final)
    # for i in range(5):
    zx.simplify.clifford_simp(pyzx_final)
    print("After clifford_simp")
    draw(pyzx_final)

    pyzx_final = bad_heuristic_are_phases_zero(pyzx_final, outer_symbol)
    print("After bad phase heuristic")
    draw(pyzx_final)
  
    zx.simplify.clifford_simp(pyzx_final)
    print("After clifford_simp")
    draw(pyzx_final)
    # zx.simplify.id_simp(pyzx_final)
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
