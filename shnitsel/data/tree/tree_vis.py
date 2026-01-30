from functools import lru_cache, partial
from importlib.resources import files
from typing import TYPE_CHECKING, Any, Callable, Hashable, Literal, Mapping
import uuid
from html import escape
from .node import TreeNode
from .datatree_level import DataTreeLevelMap
from shnitsel.bridges import default_mol
import rdkit

STATIC_FILES = (
    ("shnitsel.vis.static.html", "icons-svg-inline.html"),
    ("shnitsel.vis.static.css", "style.css"),
)

# if TYPE_CHECKING:
#     from .node import TreeNode


@lru_cache(None)
def _load_static_files():
    """Lazily load the resource files into memory the first time they are needed"""
    return [
        files(package).joinpath(resource).read_text(encoding="utf-8")
        for package, resource in STATIC_FILES
    ]


def _icon(
    icon_name: Literal[
        'database',
        'file-text2',
        'folder',
        'tree_struct',
        'tree_database',
        'molecule',
        'molecule2',
        'tree',
        'file_dark',
        'file_bright',
        'file_bright2',
        'folder_open',
        'folder_closed',
    ],
) -> str:
    # icon_name should be defined in shnitsel/vis/static/html/icon-svg-inline.html
    return f"<svg class='icon st-{icon_name}'><use xlink:href='#icon-{icon_name}'></use></svg>"


def collapsible_section(
    name, inline_details="", details="", n_items=None, enabled=True, collapsed=False
) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())

    has_items = n_items is not None and n_items
    n_items_span = "" if n_items is None else f" <span>({n_items})</span>"
    enabled = "" if enabled and has_items else "disabled"
    collapsed = "" if collapsed or not has_items else "checked"
    tip = " title='Expand/collapse section'" if enabled else ""

    return (
        f"<input id='{data_id}' class='st-section-summary-in' "
        f"type='checkbox' {enabled} {collapsed}>"
        f"<label for='{data_id}' class='st-section-summary' {tip}>"
        f"{name}:{n_items_span}</label>"
        f"<div class='st-section-inline-details'>{inline_details}</div>"
        f"<div class='st-section-details'>{details}</div>"
    )


def _mapping_section(
    mapping:Mapping[Hashable, Any],
    name:str,
    details_func:Callable,
    max_items_collapse: int = 5,
    enabled:bool=True,
    expanded:bool = False
) -> str:
    n_items = len(mapping)
    collapsed = not expanded

    inline_details = ""
    if n_items > max_items_collapse:
        inline_details = f"({max_items_collapse}/{n_items})"

    return collapsible_section(
        name,
        inline_details=inline_details,
        details=details_func(mapping),
        n_items=n_items,
        enabled=enabled,
        collapsed=collapsed,
    )


def summarize_attrs(attrs) -> str:
    attrs_dl = "".join(
        f"<dt><span>{escape(str(k))} :</span></dt><dd>{escape(str(v))}</dd>"
        for k, v in attrs.items()
    )

    return f"<dl class='st-attrs'>{attrs_dl}</dl>"


def summarize_tree_children(children: Mapping[str, TreeNode]) -> str:
    n_children = len(children)

    children_html = []
    for i, child in enumerate(children.values()):
        # if i < ceil(MAX_CHILDREN / 2) or i >= ceil(n_children - MAX_CHILDREN / 2):
        #     is_last = i == (n_children - 1)
        #     children_html.append(datatree_child_repr(child, end=is_last))
        # elif n_children > MAX_CHILDREN and i == ceil(MAX_CHILDREN / 2):
        # children_html.append("<div>...</div>")
        children_html.append(datatree_child_repr(child, end=i == (n_children - 1)))

    return "".join(
        [
            "<div style='display: inline-grid; grid-template-columns: 100%; grid-column: 1 / -1'>",
            "".join(children_html),
            "</div>",
        ]
    )


attr_section = partial(
    _mapping_section,
    name=_icon('file-text2')+"&nbsp;Attributes",
    details_func=summarize_attrs,
    max_items_collapse=10,
)

root_children_section = partial(
    _mapping_section,
    name=_icon('molecule2')+"&nbsp;Compounds",
    details_func=summarize_tree_children,
    expanded=True,
)

groups_section = partial(
    _mapping_section,
    name=_icon('tree_struct')+"&nbsp;Groups",
    details_func=summarize_tree_children,
)

leaves_section = partial(
    _mapping_section,
    name=_icon('file_dark')+"Data",
    details_func=summarize_tree_children,
)


def _obj_repr(obj, header_components, sections):
    """Return HTML repr of a shnitsel tools object.

    If CSS is not injected (untrusted notebook), fallback to the plain text repr.

    """

    header = f"<div class='st-header'>{''.join(h for h in header_components)}</div>"
    sections = "".join(f"<li class='st-section-item'>{s}</li>" for s in sections)

    icons_svg, css_style = _load_static_files()

    return (
        "<div>"
        f"{icons_svg}<style>{css_style}</style>"
        f"<pre class='st-text-repr-fallback'>{escape(repr(obj))}</pre>"
        "<div class='st-wrap' style='display:none'>"
        f"{header}"
        f"<ul class='st-sections'>{sections}</ul>"
        "</div>"
        "</div>"
    )


def datatree_child_repr(node: TreeNode, end: bool = False) -> str:
    # Wrap Tree HTML representation with a tee to the left of it.
    #
    # Enclosing HTML tag is a <div> with :code:`display: inline-grid` style.
    #
    # Turns:
    # [    title    ]
    # |   details   |
    # |_____________|
    #
    # into (A):
    # |─ [    title    ]
    # |  |   details   |
    # |  |_____________|
    #
    # or (B):
    # └─ [    title    ]
    #    |   details   |
    #    |_____________|
    end = bool(end)
    height = "100%" if end is False else "1.2em"  # height of line

    path = escape(node.path)
    sections = tree_node_sections(node, root=False)
    section_items = "".join(f"<li class='st-section-item'>{s}</li>" for s in sections)

    if node.is_level('CompoundGroup'):
        mol_data = list(node.collect_data())
        if len(mol_data)>0:
            mol_source = mol_data[0]
            try:
                mol_object = default_mol(mol_source)
                mol = rdkit.Chem.Draw.rdMolDraw2D.MolToSVG(mol_object, width=100, height=60)
                mol = ''.join(mol.split('\n')[1:])
            except: #ValueError:
                # default_mol didn't work
                mol = ''
    else:
        mol = ''

    # TODO: Can we make the group name clickable to toggle the sections below?
    # This looks like it would require the input/label pattern used above.
    html = f"""
        <div class='st-group-box'>
            <div class='st-group-box-vline' style='height: {height}'></div>
            <div class='st-group-box-hline'></div>
            <div class='st-group-box-contents'>
                <div class='st-header'>
                    <div class='st-group-name'>{path}</div>
                </div>
                {mol}
                <ul class='st-sections'>
                    {section_items}
                </ul>
            </div>
        </div>
    """
    return "".join(t.strip() for t in html.split("\n"))


def tree_repr(node: TreeNode) -> str:
    header_components = [
        f"<div class='st-obj-type'>{_icon('tree_database')}&nbsp;shnitsel.{type(node).__name__} (Level: {node._level_name or '?'}) </div>",
    ]
    if node.name is not None:
        name = escape(repr(node.name))
        header_components.append(f"<div class='st-obj-name'>{name}</div>")

    sections = tree_node_sections(node, root=True)
    return _obj_repr(node, header_components, sections)


def tree_node_sections(node: TreeNode, root: bool = False) -> list[str]:
    sections = []

    if node.children:
        children_max_items = 6
        if node._level_name == DataTreeLevelMap['root']:
            sections.append(
                root_children_section(
                    node.children, max_items_collapse=children_max_items
                )
            )
        else:
            group_children = {
                k: v
                for k, v in node.children.items()
                if v._level_name == DataTreeLevelMap['group']
            }
            data_children = {
                k: v
                for k, v in node.children.items()
                if v._level_name == DataTreeLevelMap['data']
            }
            if group_children:
                sections.append(
                    groups_section(
                        group_children, max_items_collapse=children_max_items
                    )
                )
            if data_children:
                sections.append(
                    leaves_section(data_children, max_items_collapse=children_max_items)
                )

    if node.attrs:
        sections.append(attr_section(node.attrs))

    return sections
