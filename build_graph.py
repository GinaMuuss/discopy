from discopy.drawing import draw
from discopy.quantum.zx import Z, X, Id, SWAP
from discopy.quantum.rewrite_qaoa import simplify_qaoa
import sympy
import pyzx
import itertools
import networkx
from typing import List, Tuple, Set, Dict

def get_kronecker(amount, spider):
    diagram = spider
    for i in range(1, amount):
        diagram = diagram @ spider
    return diagram


def build_diag_from_graph(num_vertices: int, edges: List[Tuple[int, int]], symbol_for_gadgets, with_input: bool = False):
    # Start with the Z spider states.
    diagram = get_kronecker(num_vertices, Z(0,1) if not with_input else Z(1,1))

    # now we built the phase gadgets one by one
    for s, t in edges:
        assert s < num_vertices, f"edge invalid {s}, {num_vertices}"
        assert t < num_vertices, f"edge invalid {t}, {num_vertices}"
        if s > t:
            s, t = t, s
        assert s < t, "edge invalid"

        print("\nAdding edge", s, t)

        # first a layer with the start index
        layer = Id(s) @ Z(1, 2) @ Id(num_vertices - s - 1)
        diagram = diagram >> layer

        # now the X spider for the phase gadget
        layer = Id(s + 1) @ X(1, 2) @ Id(num_vertices - s- 1)
        diagram = diagram >> layer

        # now add the phase spider 
        layer = Id(s + 1) @ Z(1, 0, symbol_for_gadgets) @ Id(num_vertices - s)
        diagram = diagram >> layer

        # now swap down
        would_swap_to = s + 1
        while t > would_swap_to:
            layer = Id(would_swap_to) @ SWAP @ Id(num_vertices - would_swap_to - 1)
            diagram = diagram >> layer
            would_swap_to += 1

        # now combine the phase gadget back into
        layer = Id(t) @ Z(2, 1) @ Id(num_vertices - t - 1)
        diagram = diagram >> layer

        # add identity layer
        layer = Id(num_vertices )
        diagram = diagram >> layer

    return diagram

def get_middle(num_vertices, symbol, s, t):
    assert s < num_vertices, f"edge invalid {s}, {num_vertices}"
    assert t < num_vertices, f"edge invalid {t}, {num_vertices}"
    if s > t:
        s, t = t, s
    assert s < t, "edge invalid"
    middle_part =  get_kronecker(num_vertices, X(1, 1, -symbol))
    x = Id(s) @ Z(1, 1, phase=0.5) @ Id(t - s - 1) @ Z(1, 1, phase=0.5) @ Id(num_vertices - 1 - t)
    middle_part = middle_part >> x
    middle_part = middle_part >> get_kronecker(num_vertices, X(1, 1, symbol))
    return middle_part

def get_circuit_from_graph(num_vertices: int, edges: List[Tuple[int, int]], s, t, beta, phi, with_inputs: bool = False):
    diagram = build_diag_from_graph(num_vertices, edges, -phi, with_inputs)
    m = get_middle(num_vertices, beta, s, t)
    diagram = diagram >> m
    diagram = diagram >> build_diag_from_graph(num_vertices, [(num_vertices -s -1, num_vertices -t-1) for s,t in edges], phi, with_inputs).transpose()
    return diagram

def get_neighbourhood(num_vertices: int, edges: List[Tuple[int, int]], start: int, n: int):
    queue = [(start,0)]
    visited = set([start])
    distances: Dict[int, List[int]] = {i: [] for i in range(n+1)}
    distances[0] = [start]
    output = []
    while len(queue) > 0:
        cur, distance = queue[0]
        output.append(cur)
        if distance > n:
            break
        del queue[0]
        for e in edges:
            if cur not in e:
                continue
            other = e[0] if e[0] != cur else e[1]
            if other not in visited:
                queue.append((other, distance + 1))
                visited.add(other)
                distances[distance + 1].append(other)
    print(output)
    print(distances)
    return distances[n]

def get_circuit_from_graph_interation_p(p: int, num_vertices: int, edges: List[Tuple[int, int]], s, t, betas: List, phis: List):
    assert len(betas) == p, f"length of betas should be {p} but is {len(betas)}"
    assert len(phis) == p, f"length of phis should be {p} but is {len(betas)}"
    diagram = get_circuit_from_graph(num_vertices, edges, s, t, betas[0], phis[0], with_inputs=True)
    diagram.draw()
    for i in range(1, p):
        diagram = get_kronecker(num_vertices, X(1, 1, -betas[i])) >> diagram >>  get_kronecker(num_vertices, X(1, 1, betas[i])) 

        e = []
        for j in range(num_vertices):
            e += [(j, x) for x in get_neighbourhood(num_vertices, edges, j, i + 1) if (x, j) not in e]

        diagram = build_diag_from_graph(num_vertices, e, -phis[i], with_input=True) >> diagram >> build_diag_from_graph(num_vertices, e, phis[i], with_input=True).transpose()

    return diagram

if __name__ == "__main__":
    beta = sympy.Symbol("beta", complex=True)
    phi = sympy.Symbol("phi", complex=True)
    # our favorout graph  (1,2), (1,3),(2,3), (0,1)
    n = 6
    possible_edges = list(itertools.combinations(range(n), 2))
    print(possible_edges)
    not_working = []
    #for num_edges in range(len(possible_edges) + 1):
        #for edgeset in itertools.combinations(possible_edges, num_edges):
    index_first_graph_with_n = 0
    while networkx.graph_atlas(index_first_graph_with_n).number_of_nodes() != n:
        index_first_graph_with_n += 1
    index_first_graph_with_n += 150
    index_first_graph_wiht_n_plus_one = index_first_graph_with_n
    while networkx.graph_atlas(index_first_graph_wiht_n_plus_one).number_of_nodes() != n + 1:
        index_first_graph_wiht_n_plus_one += 1
    index_first_graph_wiht_n_plus_one = index_first_graph_with_n + 1
    print("this many graphs to test ", index_first_graph_wiht_n_plus_one - index_first_graph_with_n)
    j = 0
    for i in range(index_first_graph_with_n, index_first_graph_wiht_n_plus_one):
            print("working on i", i, "that is index", i - index_first_graph_with_n )
            g = networkx.graph_atlas(i)
            if(g.number_of_nodes() < n):
                continue
            if(g.number_of_nodes() > n):
                break
            edgeset = [(s,t) for s,t  in g.edges()]
            print("edgeset", edgeset)
            for edge_to_contribute in edgeset:
                j += 1
                print(edgeset, edge_to_contribute)
                #if (1,2) not in edgeset:
                #    continue
                d = get_circuit_from_graph(n, edgeset, edge_to_contribute[0], edge_to_contribute[1], beta, phi)
                #pyzx.drawing.draw(d.to_pyzx())
                dia, scalar = simplify_qaoa(d, beta, phi)
                print()
                print(dia)
                print(scalar)
                assert len(dia.boxes) == 1
                if any(not isinstance(term, Id) for term in dia.boxes[0].terms):
                    not_working.append((edgeset, edge_to_contribute, [len(term.to_pyzx().vertices()) for term in dia.boxes[0].terms]))
    print("in total there were ", j, "tries")
    print("it did not work for ", len(not_working))
    for eset, e, count_rem in not_working:
        print(f"{eset=} {e=} {count_rem=}")
    
"""eset=[(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)] e=(0, 2) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(0, 1) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(0, 2) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(0, 3) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(1, 2) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(1, 3) count_rem=[0, 0, 0, 5]
eset=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] e=(2, 3) count_rem=[0, 0, 0, 5]
"""
