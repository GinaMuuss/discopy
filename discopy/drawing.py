# -*- coding: utf-8 -*-
""" Drawing module. """

import os
from dataclasses import dataclass
from abc import ABC, abstractmethod
from tempfile import NamedTemporaryFile, TemporaryDirectory

import networkx as nx
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch


# Mapping from attributes to mappings from box to default value.
DRAWING_ATTRIBUTES = {
    "draw_as_wires": lambda _: False,
    "draw_as_spider": lambda _: False,
    "shape": lambda _: "circle",
    "color": lambda box:
        "red" if getattr(box, "draw_as_spider", False) else "white",
    "drawing_name": lambda box: box.name,
    "tikzstyle_name": lambda box: box.name,
    "bubble_opening": lambda _: False,
    "bubble_closing": lambda _: False,
}


@dataclass
class DEFAULT:
    """ Drawing defaults. """
    aspect = "auto"
    fontsize = 12
    margins = (.05, .05)
    textpad = (.1, .1)
    color = 'white'
    use_tikzstyles = False


@dataclass
class COLORS:
    """ Drawing colours. """
    white = '#ffffff'
    red = '#e8a5a5'
    green = '#d8f8d8'
    blue = '#776ff3'
    yellow = '#f7f700'
    black = '#000000'


@dataclass
class SHAPES:
    """ Drawing shapes. """
    rectangle = 's'
    circle = 'o'


class Node:
    """ Node in a :class:`networkx.Graph`, can hold arbitrary data. """
    def __init__(self, kind, **data):
        self.kind, self.data = kind, data
        for key, value in data.items():
            setattr(self, key, value)

    def __eq__(self, other):
        return isinstance(other, Node)\
            and (self.kind, self.data) == (other.kind, other.data)

    def __repr__(self):
        return "Node({}, {})".format(repr(self.kind), ", ".join(
            "{}={}".format(key, repr(value))
            for key, value in sorted(self.data.items())))

    def __hash__(self):
        return hash(repr(self))

    __str__ = __repr__


def open_bubbles(diagram):
    """ Replace bubbles by diagrams with opening and closing boxes. """
    from discopy.monoidal import Ob, Ty, Box, Id, Functor, Bubble
    if not any(isinstance(box, Bubble) for box in diagram.boxes):
        return diagram

    class OpenBubbles(Functor):
        """ Sends a bubble to open >> Id(label) @ bubble @ Id("") >> clos. """
        def __call__(self, diagram):
            if isinstance(diagram, Bubble):
                obj = Ob(diagram.drawing_name)
                obj.draw_as_box = True
                _left, _right = Ty(obj), Ty("")
                _open = Box(
                    "_open", diagram.dom, _left @ diagram.inside.dom @ _right)
                _close = Box(
                    "_close", _left @ diagram.inside.cod @ _right, diagram.cod)
                _open.draw_as_wires, _close.draw_as_wires = True, True
                if len(diagram.dom) == len(diagram.inside.dom):
                    _open.bubble_opening = True
                if len(diagram.cod) == len(diagram.inside.cod):
                    _close.bubble_closing = True
                inside = self(diagram.inside)
                return _open >> Id(_left) @ inside @ Id(_right) >> _close
            return super().__call__(diagram)
    return OpenBubbles(lambda x: x, lambda f: f)(diagram)


def diagram_to_nx(diagram):
    """
    Builds a networkx graph, called by :meth:`Diagram.draw`.

    Parameters
    ----------
    diagram : discopy.monoidal.Diagram
        any diagram.

    Returns
    -------
    graph : networkx.Graph
        with nodes for inputs, outputs, boxes and wires.

    positions : Mapping[Node, Tuple[float, float]]
        from nodes to pairs of floats.
    """
    diagram = open_bubbles(diagram.downgrade())
    graph, pos = nx.DiGraph(), dict()

    def add_node(node, position):
        graph.add_node(node)
        pos.update({node: position})

    def add_box(scan, box, off, depth, x_pos):
        bubble = box.bubble_opening or box.bubble_closing
        node = Node("box", box=box, depth=depth)
        add_node(node, (x_pos, len(diagram) - depth - .5))
        for i, obj in enumerate(box.dom):
            wire, position = Node("dom", obj=obj, i=i, depth=depth), (
                pos[scan[off + i]][0], len(diagram) - depth - .25)
            add_node(wire, position)
            graph.add_edge(scan[off + i], wire)
            if not bubble or box.bubble_closing and i in [0, len(box.dom) - 1]:
                graph.add_edge(wire, node)
        for i, obj in enumerate(box.cod):
            wire, position = Node("cod", obj=obj, i=i, depth=depth), (
                x_pos - len(box.cod[1:]) / 2 + i, len(diagram) - depth - .75)
            add_node(wire, position)
            if not bubble or box.bubble_opening and i in [0, len(box.cod) - 1]:
                graph.add_edge(node, wire)
        if box.bubble_opening:
            for i, obj in enumerate(box.dom):
                source = Node("dom", obj=obj, i=i, depth=depth)
                target = Node("cod", obj=obj, i=i + 1, depth=depth)
                graph.add_edge(source, target)
        if box.bubble_closing:
            for i, obj in enumerate(box.cod):
                source = Node("dom", obj=obj, i=i + 1, depth=depth)
                target = Node("cod", obj=obj, i=i, depth=depth)
                graph.add_edge(source, target)
        return scan[:off]\
            + [Node("cod", obj=obj, i=i, depth=depth)
               for i, obj in enumerate(box.cod)]\
            + scan[off + len(box.dom):]

    def make_space(scan, box, off):
        if not scan:
            return 0
        half_width = len(box.cod[:-1]) / 2 + 1
        if not box.dom:
            if not off:
                x_pos = pos[scan[0]][0] - half_width
            elif off == len(scan):
                x_pos = pos[scan[-1]][0] + half_width
            else:
                right = pos[scan[off + len(box.dom)]][0]
                x_pos = (pos[scan[off - 1]][0] + right) / 2
        else:
            right = pos[scan[off + len(box.dom) - 1]][0]
            x_pos = (pos[scan[off]][0] + right) / 2
        if off and pos[scan[off - 1]][0] > x_pos - half_width:
            limit = pos[scan[off - 1]][0]
            pad = limit - x_pos + half_width
            for node, position in pos.items():
                if position[0] <= limit:
                    pos[node] = (pos[node][0] - pad, pos[node][1])
        if off + len(box.dom) < len(scan)\
                and pos[scan[off + len(box.dom)]][0] < x_pos + half_width:
            limit = pos[scan[off + len(box.dom)]][0]
            pad = x_pos + half_width - limit
            for node, position in pos.items():
                if position[0] >= limit:
                    pos[node] = (pos[node][0] + pad, pos[node][1])
        return x_pos
    scan = []
    for i, obj in enumerate(diagram.dom):
        node = Node("input", obj=obj, i=i)
        add_node(node, (i, len(diagram) or 1))
        scan.append(node)
    for depth, (box, off) in enumerate(zip(diagram.boxes, diagram.offsets)):
        x_pos = make_space(scan, box, off)
        scan = add_box(scan, box, off, depth, x_pos)
    for i, obj in enumerate(diagram.cod):
        node = Node("output", obj=obj, i=i)
        add_node(node, (pos[scan[i]][0], 0))
        graph.add_edge(scan[i], node)
    return graph, pos


class Backend(ABC):
    """ Abstract drawing backend. """
    def __init__(self):
        self.max_width = 0

    @abstractmethod
    def draw_text(self, text, i, j, **params):
        """ Draws a piece of text at a given position. """
        self.max_width = max(self.max_width, i)

    @abstractmethod
    def draw_polygon(self, *points, color=DEFAULT.color):
        """ Draws a polygon given a list of points. """
        self.max_width = max(self.max_width, max(i for i, _ in points))

    @abstractmethod
    def draw_wire(self, source, target, bend_out=False, bend_in=False):
        """ Draws a wire from source to target, possibly with a Bezier. """
        self.max_width = max(self.max_width, source[0], target[0])

    @abstractmethod
    def draw_spiders(self, graph, positions, draw_box_labels=True):
        """ Draws a list of boxes depicted as spiders. """
        spider_widths = [
            positions[n][0] for n in graph.nodes
            if n.kind == 'box' and n.box.draw_as_spider]
        if spider_widths:
            self.max_width = max(self.max_width, max(spider_widths))

    @abstractmethod
    def output(self, path=None, show=True, **params):
        """ Output the drawing. """


class TikzBackend(Backend):
    """ Tikz drawing backend. """
    def __init__(self, use_tikzstyles=None):
        self.use_tikzstyles = DEFAULT.use_tikzstyles\
            if use_tikzstyles is None else use_tikzstyles
        self.node_styles, self.edge_styles = [], []
        self.nodes, self.nodelayer, self.edgelayer = {}, [], []
        super().__init__()

    @staticmethod
    def format_color(color):
        """ Formats a color. """
        hexcode = getattr(COLORS, color)
        rgb = [
            int(hex, 16) for hex in [hexcode[1:3], hexcode[3:5], hexcode[5:]]]
        return "{{rgb,255: red,{}; green,{}; blue,{}}}".format(*rgb)

    def add_node(self, i, j, text=None, options=None):
        """ Add a node to the tikz picture, return its unique id. """
        node = len(self.nodes) + 1
        self.nodelayer.append(
            "\\node [{}] ({}) at ({}, {}) {{{}}};\n".format(
                options or "", node, i, j, text or ""))
        self.nodes.update({(i, j): node})
        return node

    def draw_text(self, text, i, j, **params):
        options = "style=none"
        if params.get("verticalalignment", "center") == "top":  # wire labels
            options += ", right"
        if 'fontsize' in params and params['fontsize'] is not None:
            options += ", scale={}".format(params['fontsize'])
        self.add_node(i, j, text, options)
        super().draw_text(text, i, j, **params)

    def draw_polygon(self, *points, color=DEFAULT.color):
        nodes = []
        for point in points:
            nodes.append(self.add_node(*point))
        nodes.append(nodes[0])
        if self.use_tikzstyles:
            style_name = "box" if color == DEFAULT.color\
                else "{}_box".format(color)
            style = "\\tikzstyle{{{}}}=[-, fill={}]\n"\
                .format(style_name, self.format_color(color))
            if style not in self.edge_styles:
                self.edge_styles.append(style)
            options = "style={}".format(style_name)
        else:
            options = "-, fill={{{}}}".format(color)
        self.edgelayer.append("\\draw [{}] {};\n".format(options, " to ".join(
            "({}.center)".format(node) for node in nodes)))
        super().draw_polygon(*points, color=color)

    def draw_wire(self, source, target, bend_out=False, bend_in=False):
        out = -90 if not bend_out or source[0] == target[0]\
            else (180 if source[0] > target[0] else 0)
        inp = 90 if not bend_in or source[0] == target[0]\
            else (180 if source[0] < target[0] else 0)
        cmd = "\\draw [in={}, out={}] ({}.center) to ({}.center);\n"
        if source not in self.nodes:
            self.add_node(*source)
        if target not in self.nodes:
            self.add_node(*target)
        self.edgelayer.append(cmd.format(
            inp, out, self.nodes[source], self.nodes[target]))
        super().draw_wire(source, target, bend_out=bend_out, bend_in=bend_in)

    def draw_spiders(self, graph, positions, draw_box_labels=True):
        spiders = [(node, node.box.color, node.box.shape)
                   for node in graph.nodes
                   if node.kind == "box" and node.box.draw_as_spider]
        for node, color, shape in spiders:
            i, j = positions[node]
            text = node.box.drawing_name if draw_box_labels else ""
            if self.use_tikzstyles:
                style = "\\tikzstyle{{{}}}=[fill={}]\n"\
                    .format(node.box.tikzstyle_name, self.format_color(color))
                if style not in self.node_styles:
                    self.node_styles.append(style)
                options = "style={}".format(node.box.tikzstyle_name)
            else:
                options = "{}, fill={}".format(shape, color)
            self.add_node(i, j, text, options)
        super().draw_spiders(graph, positions, draw_box_labels)

    def output(self, path=None, show=True, **params):
        baseline = params.get("baseline", 0)
        tikz_options = params.get("tikz_options", None)
        output_tikzstyle = self.use_tikzstyles\
            and params.get("output_tikzstyle", True)
        options = "baseline=(0.base)" if tikz_options is None\
            else "baseline=(0.base), " + tikz_options
        begin = ["\\begin{{tikzpicture}}[{}]\n".format(options)]
        nodes = ["\\begin{pgfonlayer}{nodelayer}\n",
                 "\\node (0) at (0, {}) {{}};\n".format(baseline)]\
            + self.nodelayer + ["\\end{pgfonlayer}\n"]
        edges = ["\\begin{pgfonlayer}{edgelayer}\n"] + self.edgelayer +\
                ["\\end{pgfonlayer}\n"]
        end = ["\\end{tikzpicture}\n"]
        if path is not None:
            if output_tikzstyle:
                style_path = '.'.join(path.split('.')[:-1]) + '.tikzstyles'
                with open(style_path, 'w+') as file:
                    file.writelines(["% Node styles\n"] + self.node_styles)
                    file.writelines(["% Edge styles\n"] + self.edge_styles)
            with open(path, 'w+') as file:
                file.writelines(begin + nodes + edges + end)
        elif show:  # pragma: no cover
            if output_tikzstyle:
                print(''.join(self.node_styles + self.edge_styles))
            print(''.join(begin + nodes + edges + end))


class MatBackend(Backend):
    """ Matplotlib drawing backend. """
    def __init__(self, axis=None, figsize=None):
        self.axis = axis or plt.subplots(figsize=figsize)[1]
        super().__init__()

    def draw_text(self, text, i, j, **params):
        params['fontsize'] = params.get('fontsize', None) or DEFAULT.fontsize
        self.axis.text(i, j, text, **params)
        super().draw_text(text, i, j, **params)

    def draw_polygon(self, *points, color=DEFAULT.color):
        codes = [Path.MOVETO]
        codes += len(points[1:]) * [Path.LINETO] + [Path.CLOSEPOLY]
        path = Path(points + points[:1], codes)
        self.axis.add_patch(PathPatch(path, facecolor=getattr(COLORS, color)))
        super().draw_polygon(*points, color=color)

    def draw_wire(self, source, target, bend_out=False, bend_in=False):
        mid = (target[0], source[1]) if bend_out else (source[0], target[1])
        path = Path([source, mid, target],
                    [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        self.axis.add_patch(PathPatch(path, facecolor='none'))
        super().draw_wire(source, target, bend_out=bend_out, bend_in=bend_in)

    def draw_spiders(self, graph, positions, draw_box_labels=True):
        nodes = {node for node in graph.nodes
                 if node.kind == "box" and node.box.draw_as_spider}
        shapes = {node: node.box.shape for node in nodes}
        for shape in set(shapes.values()):
            colors = {n: n.box.color for n, s in shapes.items() if s == shape}
            nodes, colors = zip(*colors.items())
            nx.draw_networkx_nodes(
                graph, positions, nodelist=nodes,
                node_color=[getattr(COLORS, color) for color in colors],
                node_shape=getattr(SHAPES, shape), ax=self.axis)
            if draw_box_labels:
                labels = {node: node.box.drawing_name for node in nodes}
                nx.draw_networkx_labels(graph, positions, labels)
        super().draw_spiders(graph, positions, draw_box_labels)

    def output(self, path=None, show=True, **params):
        xlim, ylim = params.get("xlim", None), params.get("ylim", None)
        margins = params.get("margins", DEFAULT.margins)
        aspect = params.get("aspect", DEFAULT.aspect)
        plt.margins(*margins)
        plt.subplots_adjust(
            top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        self.axis.set_aspect(aspect)
        plt.axis('off')
        if xlim is not None:
            self.axis.set_xlim(*xlim)
        if ylim is not None:
            self.axis.set_ylim(*ylim)
        if path is not None:
            plt.savefig(path)
            plt.close()
        if show:
            plt.show()


def draw(diagram, **params):
    """
    Draws a diagram using networkx and matplotlib.

    Parameters
    ----------
    draw_as_nodes : bool, optional
        Whether to draw boxes as nodes, default is :code:`False`.
    color : string, optional
        Color of the box or node, default is white (:code:`'#ffffff'`) for
        boxes and red (:code:`'#ff0000'`) for nodes.
    textpad : pair of floats, optional
        Padding between text and wires, default is :code:`(0.1, 0.1)`.
    draw_type_labels : bool, optional
        Whether to draw type labels, default is :code:`False`.
    draw_box_labels : bool, optional
        Whether to draw box labels, default is :code:`True`.
    aspect : string, optional
        Aspect ratio, one of :code:`['auto', 'equal']`.
    margins : tuple, optional
        Margins, default is :code:`(0.05, 0.05)`.
    fontsize : int, optional
        Font size for the boxes, default is :code:`12`.
    fontsize_types : int, optional
        Font size for the types, default is :code:`12`.
    figsize : tuple, optional
        Figure size.
    path : str, optional
        Where to save the image, if `None` we call :code:`plt.show()`.
    to_tikz : bool, optional
        Whether to output tikz code instead of matplotlib.
    asymmetry : float, optional
        Make a box and its dagger mirror images, default is
        :code:`.25 * any(box.is_dagger for box in diagram.boxes)`.
    """
    def draw_box(backend, positions, box, depth):
        asymmetry = params.get(
            'asymmetry', .25 * any(box.is_dagger for box in diagram.boxes))
        node = Node("box", box=box, depth=depth)
        if not box.dom and not box.cod:
            left, right = positions[node][0], positions[node][0]
        elif not box.dom:
            left, right = (
                positions[Node("cod", obj=box.cod[i], i=i, depth=depth)][0]
                for i in [0, len(box.cod) - 1])
        elif not box.cod:
            left, right = (
                positions[Node("dom", obj=box.dom[i], i=i, depth=depth)][0]
                for i in [0, len(box.dom) - 1])
        else:
            top_left, top_right = (
                positions[Node("dom", obj=box.dom[i], i=i, depth=depth)][0]
                for i in [0, len(box.dom) - 1])
            bottom_left, bottom_right = (
                positions[Node("cod", obj=box.cod[i], i=i, depth=depth)][0]
                for i in [0, len(box.cod) - 1])
            left = min(top_left, bottom_left)
            right = max(top_right, bottom_right)
        height = positions[node][1] - .25
        left, right = left - .25, right + .25
        backend.draw_polygon(
            (left, height),
            (right + (asymmetry if box.is_dagger else 0), height),
            (right + (0 if box.is_dagger else asymmetry), height + .5),
            (left, height + .5),
            color=params.get('color', DEFAULT.color))
        if params.get('draw_box_labels', True):
            backend.draw_text(box.drawing_name, *positions[node],
                              ha='center', va='center',
                              fontsize=params.get('fontsize', None))
        return backend

    def draw_wires(backend, graph, positions):
        for source, target in graph.edges():
            def inside_a_box(node):
                return node.kind == "box"\
                    and not node.box.draw_as_wires\
                    and not node.box.draw_as_spider
            if inside_a_box(source) or inside_a_box(target):
                continue  # no need to draw wires inside a box
            backend.draw_wire(
                positions[source], positions[target],
                bend_out=source.kind == "box", bend_in=target.kind == "box")
            if source.kind in ["input", "cod"]\
                    and (params.get('draw_type_labels', True)
                         or getattr(source.obj, "draw_as_box", False)
                         and params.get('draw_box_labels', True)):
                i, j = positions[source]
                pad_i, pad_j = params.get('textpad', DEFAULT.textpad)
                pad_j = 0 if source.kind == "input" else pad_j
                backend.draw_text(
                    str(source.obj), i + pad_i, j - pad_j,
                    fontsize=params.get('fontsize_types',
                                        params.get('fontsize', None)),
                    verticalalignment='top')
        return backend

    def scale_and_pad(graph, pos, scale, pad):
        widths, heights = zip(*pos.values())
        min_width, min_height = min(widths), min(heights)
        pos = {n: ((x - min_width) * scale[0] + pad[0],
                   (y - min_height) * scale[1] + pad[1])
               for n, (x, y) in pos.items()}
        for box_node in graph.nodes:
            if box_node.kind == "box":
                for i, obj in enumerate(box_node.box.dom):
                    node = Node("dom", obj=obj, i=i, depth=box_node.depth)
                    pos[node] = (
                        pos[node][0], pos[node][1] - .25 * (scale[1] - 1))
                for i, obj in enumerate(box_node.box.cod):
                    node = Node("cod", obj=obj, i=i, depth=box_node.depth)
                    pos[node] = (
                        pos[node][0], pos[node][1] + .25 * (scale[1] - 1))
        return pos

    scale, pad = params.get('scale', (1, 1)), params.get('pad', (0, 0))
    graph, positions = diagram_to_nx(diagram)
    positions = scale_and_pad(graph, positions, scale, pad)
    backend = params['backend'] if 'backend' in params else\
        TikzBackend(use_tikzstyles=params.get('use_tikzstyles', None))\
        if params.get('to_tikz', False)\
        else MatBackend(figsize=params.get('figsize', None))

    backend = draw_wires(backend, graph, positions)
    backend.draw_spiders(
        graph, positions,
        draw_box_labels=params.get('draw_box_labels', True))
    box_nodes = [node for node in graph.nodes if node.kind == "box"]
    for node in box_nodes:
        if node.box.draw_as_spider or node.box.draw_as_wires:
            continue
        backend = draw_box(backend, positions, node.box, node.depth)
    return backend.output(
        path=params.get('path', None),
        baseline=len(box_nodes) / 2 or .5,
        tikz_options=params.get('tikz_options', None),
        show=params.get('show', True),
        margins=params.get('margins', DEFAULT.margins),
        aspect=params.get('aspect', DEFAULT.aspect))


def to_gif(diagram, *diagrams, **params):  # pragma: no cover
    """
    Builds a gif with the normalisation steps.

    Parameters
    ----------
    diagrams : :class:`Diagram`, optional
        Sequence of diagrams to draw.
    path : str
        Where to save the image, if :code:`None` a gif gets created.
    timestep : int, optional
        Time step in milliseconds, default is :code:`500`.
    loop : bool, optional
        Whether to loop, default is :code:`False`
    params : any, optional
        Passed to :meth:`Diagram.draw`.
    """
    path = params.get("path", None)
    timestep = params.get("timestep", 500)
    loop = params.get("loop", False)
    steps, frames = (diagram, ) + diagrams, []
    path = path or os.path.basename(NamedTemporaryFile(
        suffix='.gif', prefix='tmp_', dir='.').name)
    with TemporaryDirectory() as directory:
        for i, _diagram in enumerate(steps):
            tmp_path = os.path.join(directory, '{}.png'.format(i))
            _diagram.draw(path=tmp_path, **params)
            frames.append(Image.open(tmp_path))
        if loop:
            frames = frames + frames[::-1]
        frames[0].save(path, format='GIF', append_images=frames[1:],
                       save_all=True, duration=timestep,
                       **{'loop': 0} if loop else {})
        try:
            from IPython.display import HTML
            return HTML('<img src="{}">'.format(path))
        except ImportError:
            return '<img src="{}">'.format(path)


def pregroup_draw(words, cups, **params):
    """
    Draws pregroup words and cups.
    """
    textpad = params.get('textpad', (.1, .2))
    textpad_words = params.get('textpad_words', (0, .1))
    space = params.get('space', .5)
    width = params.get('width', 2.)
    fontsize = params.get('fontsize', None)

    backend = TikzBackend(use_tikzstyles=params.get('use_tikzstyles', None))\
        if params.get('to_tikz', False)\
        else MatBackend(figsize=params.get('figsize', None))

    def draw_triangles(words):
        scan = []
        for i, word in enumerate(words.boxes):
            for j, _ in enumerate(word.cod):
                x_wire = (space + width) * i\
                    + (width / (len(word.cod) + 1)) * (j + 1)
                scan.append(x_wire)
                if params.get('draw_type_labels', True):
                    backend.draw_text(
                        str(word.cod[j]), x_wire + textpad[0], -textpad[1],
                        fontsize=params.get('fontsize_types', fontsize))
            backend.draw_polygon(
                ((space + width) * i, 0),
                ((space + width) * i + width, 0),
                ((space + width) * i + width / 2, 1),
                color=DEFAULT.color)
            backend.draw_text(
                str(word), (space + width) * i + width / 2 + textpad_words[0],
                textpad_words[1], ha='center', fontsize=fontsize)
        return scan

    def draw_cups_and_wires(cups, scan):
        for j, off in [(j, off)
                       for j, s in enumerate(cups) for off in s.offsets]:
            middle = (scan[off] + scan[off + 1]) / 2
            backend.draw_wire((scan[off], 0), (middle, - j - 1), bend_in=True)
            backend.draw_wire(
                (scan[off + 1], 0), (middle, - j - 1), bend_in=True)
            scan = scan[:off] + scan[off + 2:]
        for i, _ in enumerate(cups[-1].cod if cups else words.cod):
            label = str(cups[-1].cod[i]) if cups else ""
            backend.draw_wire((scan[i], 0), (scan[i], - (len(cups) or 1) - 1))
            if params.get('draw_type_labels', True):
                backend.draw_text(
                    label, scan[i] + textpad[0], - (len(cups) or 1) - space,
                    fontsize=params.get('fontsize_types', fontsize))

    scan = draw_triangles(words.normal_form())
    draw_cups_and_wires(cups, scan)
    backend.output(
        params.get('path', None),
        tikz_options=params.get('tikz_options', None),
        xlim=(0, (space + width) * len(words.boxes) - space),
        ylim=(- len(cups) - space, 1),
        margins=params.get('margins', DEFAULT.margins),
        aspect=params.get('aspect', DEFAULT.aspect))


def equation(*diagrams, path=None, symbol="=", space=1, **params):
    """ Draws an equation with multiple diagrams. """
    def height(diagram):
        if hasattr(diagram, "terms"):  # i.e. if isinstance(diagram, Sum)
            return max(height(d) for d in diagram.terms)
        if hasattr(diagram, "inside"):  # i.e. if isinstance(diagram, Bubble)
            return height(diagram.inside) + 2
        if len(diagram) > 1:
            return sum(height(d) for d in diagram.boxes)
        return 1

    pad, max_height = params.get('pad', (0, 0)), max(map(height, diagrams))
    scale_x, scale_y = params.get('scale', (1, 1))
    backend = params['backend'] if 'backend' in params\
        else TikzBackend(use_tikzstyles=params.get('use_tikzstyles', None))\
        if params.get('to_tikz', False)\
        else MatBackend(figsize=params.get('figsize', None))

    for i, diagram in enumerate(diagrams):
        scale = (scale_x, scale_y * max_height / height(diagram))
        diagram.draw(**dict(
            params, show=False, path=None,
            backend=backend, scale=scale, pad=pad))
        pad = (backend.max_width + space, 0)
        if i < len(diagrams) - 1:
            backend.draw_text(symbol, pad[0], scale_y * max_height / 2)
            pad = (pad[0] + space, pad[1])

    return backend.output(
        path=path,
        baseline=max_height / 2,
        tikz_options=params.get('tikz_options', None),
        show=params.get("show", True),
        margins=params.get('margins', DEFAULT.margins),
        aspect=params.get('aspect', DEFAULT.aspect))


class Equation:
    """
    An equation is a list of diagrams with a dedicated draw method.

    Example
    -------
    >>> from discopy.tensor import Spider, Swap, Dim, Id
    >>> dim = Dim(2)
    >>> mu, eta = Spider(2, 1, dim), Spider(0, 1, dim)
    >>> delta, upsilon = Spider(1, 2, dim), Spider(1, 0, dim)
    >>> special = Equation(mu >> delta, Id(dim))
    >>> special  # doctest: +ELLIPSIS
    Equation(Diagram(...), Id(Dim(2)))
    >>> frobenius = Equation(
    ...     delta @ Id(dim) >> Id(dim) @ mu,
    ...     mu >> delta,
    ...     Id(dim) @ delta >> mu @ Id(dim))
    >>> print(frobenius)  # doctest: +ELLIPSIS
    Spider... @ Spider... = Spider... >> Spider... = Id... @ Spider...
    >>> equation(special, frobenius, symbol=', ',
    ...          aspect='equal', draw_type_labels=False, figsize=(8, 2),
    ...          path='docs/_static/imgs/drawing/frobenius-axioms.png')

    .. image:: ../_static/imgs/drawing/frobenius-axioms.png
        :align: center
    """
    def __init__(self, *terms, symbol='='):
        self.terms, self.symbol = terms, symbol

    def __repr__(self):
        return "Equation({})".format(', '.join(map(repr, self.terms)))

    def __str__(self):
        return " {} ".format(self.symbol).join(map(str, self.terms))

    def draw(self, **params):
        """ Drawing an equation. """
        return equation(*self.terms, **dict(params, symbol=self.symbol))
