import discopy.cat as cat
import discopy.monoidal as monoidal

def distribute_composition_cat(
    arrow: cat.Arrow,
    index_sum: int,
    distribute_up_to: int,
    arrow_construction=lambda dom, cod, boxes: cat.Arrow(dom, cod, boxes),
) -> cat.Arrow:
    """
    Distributes the sum at index_sum to the left with

    Parameters
    ----------
    arrow: Arrow
        The arrow to perform the operation on
    index_sum: int
        The index of the sum to distribute for
    distribute_up_to: int
        index of the last term to distribute
    arrow_construction: callable
        callable returning the arrow equivalent in the relevant category

    Raises
    ------
      IndexError: The index of sum does not correspond to a box
                  or distribute_up_to is negative or to large
      TypeError: The box at index_sum does not have type sum
    """
    if len(arrow.boxes) <= index_sum or index_sum < 0:
        raise IndexError("index to large or negative, no such box")
    if (
        len(arrow.boxes) <= distribute_up_to
        or distribute_up_to < 0
        or distribute_up_to == index_sum
    ):
        raise IndexError(
            "distribute_up_to to large, negative or \
equal to index_sum, no such box"
        )
    if not isinstance(arrow.boxes[index_sum], cat.LocalSum):
        raise TypeError("box at index %d is not a LocalSum", index_sum)

    if index_sum < distribute_up_to:
        # distribute to the right
        unit = arrow.boxes[index_sum].__class__(
            [], arrow.boxes[index_sum].dom, arrow.boxes[distribute_up_to].cod
        )
        terms = [
            g.then(
                arrow_construction(
                    arrow.boxes[index_sum + 1].dom,
                    arrow.boxes[distribute_up_to].cod,
                    arrow.boxes[index_sum + 1: distribute_up_to + 1],
                )
            )
            for g in arrow.boxes[index_sum].terms
        ]
        term = arrow.boxes[index_sum].upgrade(cat.LocalSum(terms))
        new_boxes = (
            arrow.boxes[:index_sum] + [term]
            + arrow.boxes[distribute_up_to + 1:]
        )
        return arrow_construction(
            new_boxes[0].dom,
            new_boxes[-1].cod,
            new_boxes,
        )
    else:
        # distribute to the left
        unit = arrow.boxes[index_sum].__class__(
            [], arrow.boxes[distribute_up_to].dom, arrow.boxes[index_sum].cod
        )
        terms = [
            arrow_construction(
                arrow.boxes[distribute_up_to].dom,
                arrow.boxes[index_sum - 1].cod,
                arrow.boxes[distribute_up_to:index_sum],
            ).then(g)
            for g in arrow.boxes[index_sum].terms
        ]
        term = arrow.boxes[index_sum].upgrade(cat.LocalSum(terms))
        new_boxes = (
            arrow.boxes[:distribute_up_to]
            + [term] + arrow.boxes[index_sum + 1:]
        )
        return arrow_construction(
            new_boxes[0].dom,
            new_boxes[-1].cod,
            new_boxes,
        )


def distribute_tensor(
    diagram: monoidal.Diagram, index_sum: int, index_of_partner: int
) -> monoidal.Diagram:
    """
    Distributes the sum at index_sum to the left with

    Parameters
    ----------
    diagram: Diagram
        The arrow to perform the operation on
    index_sum: int
        The index of the sum to distribute for
    index_of_partner: int
        index of the term to distribute

    Raises
    ------
      IndexError: The index of sum does not correspond to a box or index_of_partner is negative or to large
      TypeError: The box at index_sum does not have type sum
    """
    if len(diagram.boxes) <= index_sum or index_sum < 0:
        raise IndexError("invalid index_sum")
    if len(diagram.boxes) <= index_of_partner or index_of_partner < 0:
        raise IndexError("invalid index_of_partner")
    if index_of_partner == index_sum:
        raise IndexError("cannot distribte over itself")
    if not isinstance(diagram.boxes[index_sum], monoidal.LocalSum):
        raise TypeError("box at index_sum not LocalSum")
    if abs(index_sum - index_of_partner) != 1:
        raise IndexError(
            "at the moment distributing is only supported when the boxes are next to one another"
        )

    layers_index = list(
        range(
            diagram.offsets[index_sum],
            diagram.offsets[index_sum] + len(diagram.boxes[index_sum].dom.objects),
        )
    )
    layers_partner = list(
        range(
            diagram.offsets[index_of_partner],
            diagram.offsets[index_of_partner]
            + len(diagram.boxes[index_of_partner].dom.objects),
        )
    )

    if index_sum < index_of_partner:
        if min(layers_index) < min(layers_partner) and max(layers_index) >= min(
            layers_partner
        ):
            raise IndexError("the layers overlap! we cannot distribute")

        if min(layers_index) > min(layers_partner) and max(layers_index) <= min(
            layers_partner
        ):
            raise IndexError("the layers overlap! we cannot distribute")

    if layers_index[0] < layers_partner[0]:
        # tensor from the right
        new_terms = (
            diagram.boxes[: min(index_sum, index_of_partner)]
            + [
                monoidal.LocalSum(
                    [
                        f.tensor(diagram.boxes[index_of_partner])
                        for f in diagram.boxes[index_sum].terms
                    ]
                )
            ]
            + diagram.boxes[max(index_sum, index_of_partner) +1:]
        )
    else:
        # tensor from the left
        new_terms = (
            diagram.boxes[: min(index_sum, index_of_partner)]
            + [
                monoidal.LocalSum(
                    [
                        diagram.boxes[index_of_partner].tensor(f)
                        for f in diagram.boxes[index_sum].terms
                    ]
                )
            ]
            + diagram.boxes[max(index_sum, index_of_partner)+1 :]
        )
    new_offsets = (
        diagram.offsets[: min(index_sum, index_of_partner)]
        + [min(layers_index[0], layers_partner[0])]
        + diagram.offsets[max(index_sum, index_of_partner) +1 :]
    )
    return monoidal.Diagram(diagram.dom, diagram.cod, new_terms, new_offsets)



