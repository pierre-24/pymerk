from docutils import nodes
from docutils.nodes import section
from sphinx.util.docutils import SphinxDirective, SphinxRole

class PyMERKKeywordDef(SphinxDirective):
    """Directive for a given keyword for pymerk TOML input.

    Inspired by https://sphinx-hxl.readthedocs.io/en/stable/extdev/tutorial.html
    """
    has_content = True
    required_arguments = 1
    option_spec = {
        'type': str,
        'default': str,
        'section': str,
    }

    def run(self):
        if 'type' not in self.options:
            return [nodes.literal(text='missing `type` to .. pymkwdef')]

        if 'section' not in self.options:
            return [nodes.literal(text='missing `section` to .. pymkwdef')]

        keyword_name = '{}.{}'.format(self.options['section'], self.arguments[0])
        target_id = 'kw-{}-{}'.format(self.options['section'].lower(), self.arguments[0].lower().replace('_', ''))

        # keep the list of keywords
        if not hasattr(self.env, 'pymerk_keywords'):
            self.env.pymerk_keywords = {}

        self.env.pymerk_keywords[keyword_name] = {
            'docname': self.env.current_document.docname,
            'target': target_id,
            'section': self.options['section'],
            'name': self.arguments[0],
        }

        # make nodes
        item_node = nodes.definition_list_item()

        reference_node = nodes.reference('', '', refuri='#{}'.format(target_id))
        reference_node += nodes.literal(text='[{}]'.format(self.options['section']))
        reference_node += nodes.Text(' ')
        reference_node += nodes.strong(text=self.arguments[0], classes=['pymkw'])

        target_node = nodes.target('', '', ids=[target_id])
        target_node += reference_node
        target_node += nodes.Text(' (')
        target_node += nodes.literal(text=self.options.get('type', 'None'))

        if 'default' in self.options:
            target_node += nodes.Text(', default: ')
            target_node += nodes.literal(text=self.options.get('default'))

        target_node += nodes.Text(')')

        term_node = nodes.term()
        term_node += target_node
        item_node += term_node

        definition_node = nodes.definition()
        definition_node += self.parse_content_to_nodes()

        item_node += definition_node

        return [item_node]


class PstliteKeywordNode(nodes.General, nodes.Element):
    """Dummy node, to be ultimately replaced by a reference"""
    def __init__(self, keyword_name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.keyword_name = keyword_name


class PyMERKKeywordRole(SphinxRole):
    """Role for keyword use"""
    def run(self):
        node = PstliteKeywordNode(self.text)
        return [node], []


def process_keyword_nodes(app, doctree, fromdocname):
    """Replace `PstliteKeywordNode` by a `nodes.reference`.

    Inspired by https://sphinx-hxl.readthedocs.io/en/stable/extdev/tutorial.html
    """

    keywords = {}
    if hasattr(app.env, 'pymerk_keywords'):
        keywords = app.env.pymerk_keywords

    for node in doctree.findall(PstliteKeywordNode):
        if node.keyword_name in keywords:
            kw = keywords[node.keyword_name]

            uri = app.builder.get_relative_uri(fromdocname, kw['docname'])
            reference_node = nodes.reference('', '', refuri='{}#{}'.format(uri, kw['target']))

            reference_node += nodes.literal(text='[{}]'.format(kw['section']))
            reference_node += nodes.Text(' ')
            reference_node += nodes.strong(text=kw['name'], classes=['pymkw'])

            node.replace_self(reference_node)
        else:
            node.replace_self(nodes.literal(text='unknown keyword `{}`'.format(node.keyword_name)))


def setup(app):
    app.add_directive('pymkwdef', PyMERKKeywordDef)
    app.add_role('pymkw', PyMERKKeywordRole())

    app.connect('doctree-resolved', process_keyword_nodes)

    return {'version': '0.1'}