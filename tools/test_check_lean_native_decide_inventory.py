from tools.check_lean_native_decide_inventory import NATIVE_TOKEN_RE, code_view


def test_comments_keep_token_boundaries_and_source_locations():
    text = "by/- annotation /- nested -/\n -/native_decide\nnative_/-separator-/decide\n"
    visible = code_view(text)
    matches = list(NATIVE_TOKEN_RE.finditer(visible))
    assert [match.start() for match in matches] == [text.index("native_decide")]
    assert len(text) == len(visible)
    assert text.splitlines()[1].index("native_decide") == visible.splitlines()[1].index("native_decide")


def test_raw_and_quoted_strings_and_quote_character_are_not_tactics():
    text = """def source := r##"an internal "native_decide" and "#
still raw"##
def quotation := '"'
def message := "native_decide"
-- native_decide
example : True := by native_decide
"""
    assert [match.start() for match in NATIVE_TOKEN_RE.finditer(code_view(text))] == [
        text.rindex("native_decide")
    ]


def test_literal_contents_cannot_supply_a_declaration_or_axiom_receipt():
    from tools.check_lean_docstring_style import scan_module
    from tools.check_lean_theorem_count import DECL_RE, strip_comments

    text = '''def source := r#"
theorem phantom : True := trivial
#print axioms actual
"#
theorem actual : True := trivial
'''
    assert len(DECL_RE.findall(strip_comments(text))) == 1
    scan = scan_module(text)
    assert scan.declared == ["actual"]
    assert scan.entries == []
