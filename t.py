from discopy.drawing import draw
from discopy.quantum.zx import Z, X, Id, SWAP
from discopy.quantum.rewrite_qaoa import simplify_qaoa
from sympy.abc import phi, beta
import pyzx


def get_beginning(symbol):
    diagram = Z(0, 1) @ Z(0, 1) @ Z(0, 1) @ Z(0, 1)
    diagram = diagram >> Z(1, 2) @ Id(1) @ Z(1, 3) @ Id(1)
    diagram = diagram >> Id(1) @ X(1, 2) @ Id(1) @ X(1, 2) @ Id(1) @ X(1, 2) @ Id(1)
    diagram = diagram >> Id(7) @ Z(1, 0, symbol) @ Z(2, 2)
    diagram = diagram >> Id(7) @ X(1, 2) @ Id(1)
    diagram = diagram >> Id(1) @ Z(1, 0, symbol) @ Id(3) @ Z(1, 0, symbol) @ SWAP @ Z(
        1, 0, symbol
    ) @ Id(1)
    return diagram >> Id(1) @ Z(4, 1) @ Id(2) >> Z(1, 1) @ Z(1, 1) @ Z(1, 1) @ Z(1, 1)


def get_middle(symbol):
    middle_part = (
        X(1, 1, -symbol) @ X(1, 1, -symbol) @ X(1, 1, -symbol) @ X(1, 1, -symbol)
    )
    middle_part = middle_part >> Id(1) @ Z(1, 1, phase=0.5) @ Z(1, 1, phase=0.5) @ Id(1)
    middle_part = middle_part >> X(1, 1, symbol) @ X(1, 1, symbol) @ X(
        1, 1, symbol
    ) @ X(1, 1, symbol)
    return middle_part


def get_graph(inner_symbol, outer_symbol):
    stageone = get_beginning(-outer_symbol)
    middle_part = get_middle(inner_symbol)
    diagram = stageone >> middle_part
    final_diag = (
        diagram
        >> SWAP @ SWAP
        >> Id(1) @ SWAP @ Id(1)
        >> SWAP @ SWAP
        >> Id(1) @ SWAP @ Id(1)
        >> get_beginning(outer_symbol).transpose()
    )
    return final_diag

def main():
    a = get_graph(beta, phi)
    #b = a.to_pyzx()
    #c = b.to_tikz()
    #with open("a", "w") as f:
    #    f.write(c)
    #print(c)
    
    b = simplify_qaoa(a, beta, phi)
    b.draw()

if __name__ == "__main__":
    main()
