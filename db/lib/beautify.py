"""Converts the different BOINC result schemas to fancy TeX documents."""
import json
import os
from sympy import symbols, latex
from subprocess import DEVNULL, CalledProcessError, check_call
import PIL
from typing import Union, List
from urllib.request import Request, urlopen


TEMPLATES = {
    "conjecture": r"""\documentclass[preview,border=2pt]{standalone}
\usepackage{amsmath}
\usepackage{graphicx}
\begin{document}
\begin{equation*}
\resizebox{\textwidth}{!}
{%
    $ LHSEquation=\displaystyle\lim _{n\to\infty} RHSEquation $%
}
\end{equation*}
\end{document}""",
    "unknown-lhs": r"""\documentclass[preview,border=2pt]{standalone}
\usepackage{amsmath}
\usepackage{graphicx}
\begin{document}
\begin{equation*}
\resizebox{\textwidth}{!}
{%
    $ LHSEquation = RHSEquation $%
}
\end{equation*}
\end{document}""",
}

n = symbols("n")


def format_tex(
    template="conjecture",
    an_equation="",
    bn_equation="",
    lhs_equation="",
    rhs_equation="",
):
    return (
        TEMPLATES[template]
        .replace("AnEquation", an_equation)
        .replace("BnEquation", bn_equation)
        .replace("LHSEquation", lhs_equation)
        .replace("RHSEquation", rhs_equation)
    )


def generate_rhs(
    an_equation: str, bn_equation: str, steps=3, show_substitution=False
) -> str:
    def replace_n(equation: str, n):
        return (
            equation.replace("0n", rf"0\cdot {n}")
            .replace("1n", rf"1\cdot {n}")
            .replace("2n", rf"2\cdot {n}")
            .replace("3n", rf"3\cdot {n}")
            .replace("4n", rf"4\cdot {n}")
            .replace("5n", rf"5\cdot {n}")
            .replace("6n", rf"6\cdot {n}")
            .replace("7n", rf"7\cdot {n}")
            .replace("8n", rf"8\cdot {n}")
            .replace("9n", rf"9\cdot {n}")
            .replace("0 n", rf"0\cdot {n}")
            .replace("1 n", rf"1\cdot {n}")
            .replace("2 n", rf"2\cdot {n}")
            .replace("3 n", rf"3\cdot {n}")
            .replace("4 n", rf"4\cdot {n}")
            .replace("5 n", rf"5\cdot {n}")
            .replace("6 n", rf"6\cdot {n}")
            .replace("7 n", rf"7\cdot {n}")
            .replace("8 n", rf"8\cdot {n}")
            .replace("9 n", rf"9\cdot {n}")
            .replace("n", str(n))
        )

    def tex_to_expr(tex: str):
        return (
            tex.replace("0n", "0*n")
            .replace("1n", "1*n")
            .replace("2n", "2*n")
            .replace("3n", "3*n")
            .replace("4n", "4*n")
            .replace("5n", "5*n")
            .replace("6n", "6*n")
            .replace("8n", "8*n")
            .replace("9n", "9*n")
            .replace("0 n", "0*n")
            .replace("1 n", "1*n")
            .replace("2 n", "2*n")
            .replace("3 n", "3*n")
            .replace("4 n", "4*n")
            .replace("5 n", "5*n")
            .replace("6 n", "6*n")
            .replace("7 n", "7*n")
            .replace("8 n", "8*n")
            .replace("9 n", "9*n")
            .replace("^", "**")
            .replace(r"\cdot", "*")
            .replace("{", "(")
            .replace("}", ")")
        )

    def an(n):
        if show_substitution:
            return replace_n(an_equation, n)
        return eval(tex_to_expr(an_equation))

    def bn(n):
        if show_substitution:
            return replace_n(bn_equation, n)
        return eval(tex_to_expr(bn_equation))

    rhs_equation = (
        r"\cdots"
        if show_substitution
        else (rf"\ddots + \cfrac{{{bn_equation}}}{{{an_equation}}}")
    )
    for i in range(steps, 0, -1):
        rhs_equation = rf"{an(i-1)}+\cfrac{{{bn(i)}}}{{{rhs_equation}}}"

    return rhs_equation


def coefficient_to_tex(coefficient: int, term: str, add_dot_and_parantheses=True):
    if coefficient == 0:
        return ""
    if coefficient == 1 and term != "":
        return term
    if coefficient == -1 and term != "":
        return "-" + term
    if add_dot_and_parantheses:
        return rf"{coefficient} \cdot ({term})"
    return f"{coefficient} {term}"


def plus_coefficient_to_tex(coefficient: int, term: str, add_dot_and_parantheses=True):
    return ("+" if coefficient > 0 else "") + coefficient_to_tex(
        coefficient, term, add_dot_and_parantheses
    )


def create_consts_sum_tex(coefficients: List[int], consts: List[str]) -> str:
    assert len(coefficients) == len(consts)
    result = coefficient_to_tex(
        coefficients[0], consts[0], add_dot_and_parantheses=False
    )
    for i in range(1, len(coefficients)):
        result += plus_coefficient_to_tex(
            coefficients[i], consts[i], add_dot_and_parantheses=False
        )
    return result


def fraction(numerator: str, denominator: str) -> str:
    numerator = numerator.removeprefix("+")
    denominator = denominator.removeprefix("+")
    return rf"\cfrac{{{numerator}}}{{{denominator}}}"


def handle_general(result_data):
    # The two first lists are in the same format so they should be treated the same
    an_bn_equations = []
    for i in range(2):
        polynomial = 0
        for index, coefficient in enumerate(result_data[i][::-1]):
            polynomial += coefficient * (n**index)
        an_bn_equations.append(latex(polynomial))
    an_equation, bn_equation = an_bn_equations

    consts = ["", r"\zeta (3)", r"\zeta (2)"]
    lhs_equation = str(round(float(result_data[2]), 10))
    if result_data[3] is not None and len(result_data[3]) == 3:
        lhs_numerator = create_consts_sum_tex(result_data[3], consts)
        lhs_denominator = create_consts_sum_tex(result_data[4], consts)
        lhs_equation = fraction(lhs_numerator, lhs_denominator)

    return (
        "conjecture" if lhs_equation != "" else "unknown-lhs",
        an_equation,
        bn_equation,
        lhs_equation,
    )


def handle_zeta5(result_data):
    an_equation = (
        coefficient_to_tex(result_data[0][0], "n^5 + (n + 1)^5")
        + plus_coefficient_to_tex(result_data[0][1], "n^3 + (n + 1)^3")
        + plus_coefficient_to_tex(result_data[0][2], "2n + 1")
    )
    bn_equation = coefficient_to_tex(
        -(result_data[1][0] ** 2), "n^{10}", add_dot_and_parantheses=False
    )

    consts = ["", r"\zeta (3)", r"\zeta (5)"]

    lhs_equation = str(round(float(result_data[2]), 10))
    if result_data[3] is not None and len(result_data[3]) == 3:
        lhs_numerator = create_consts_sum_tex(result_data[3], consts)
        lhs_denominator = create_consts_sum_tex(result_data[4], consts)
        lhs_equation = fraction(lhs_numerator, lhs_denominator)

    return (
        "conjecture" if lhs_equation != "" else "unknown-lhs",
        an_equation,
        bn_equation,
        lhs_equation,
    )


HANDLERS = {"general": handle_general, "zeta5": handle_zeta5}

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def generate_tex(result_filename: str, result: str, schema: Union[str, None] = None):
    result_data = json.loads(result)[0]

    result_type = schema if schema is not None else os.path.basename(result_filename).split("_")[1]
    if result_type in HANDLERS:
        template, an_equation, bn_equation, lhs_equation = HANDLERS[result_type](
            result_data
        )
        an_equation = remove_prefix(an_equation, "+")
        bn_equation = remove_prefix(bn_equation, "+")
        rhs_equation = generate_rhs(an_equation, bn_equation)
        return (
            format_tex(template, an_equation, bn_equation, lhs_equation, rhs_equation),
            template,
        )
    else:
        print(
            f"Unsupported result type '{result_type}'\nThe supported types are: {', '.join(HANDLERS.keys())}"
        )

def execute_silently(command: str, ignore_codes: List[int] = []) -> int:
    try:
        check_call(command.split(" "), stdout=DEVNULL, stderr=DEVNULL, stdin=DEVNULL)
    except CalledProcessError as e:
        if e.returncode in ignore_codes:
            return
        print(f"Failure running '{command}': {e}")
        exit(-1)


def render_preview(
    filename: str,
    margins=2,
    density=2000,
    scale=1,
    transparent=True,
    max_aspect_ratio: Union[float, None] = None,
    pdf_filename="preview",
    cropped_pdf_filename="preview.pdf",
    preview_filename="preview.png",
) -> None:
    execute_silently(f"pdflatex -jobname={pdf_filename} {filename}")
    execute_silently(
        f"pdfcrop --margins {margins} {pdf_filename} {cropped_pdf_filename}"
    )
    transparency_option = "" if transparent else "-alpha off"
    execute_silently(
        f"convert {transparency_option} -density {density} {cropped_pdf_filename} -scale {scale * 100}% {preview_filename}",
        ignore_codes=[1],
    )

    if max_aspect_ratio:
        height = None
        width = None
        with PIL.Image.open(preview_filename) as preview_image:
            height = preview_image.height
            width = preview_image.width
        if width / height > max_aspect_ratio:
            execute_silently(
                f"convert {preview_filename} -background white -thumbnail {width}x{int(width / max_aspect_ratio)}> -gravity center -extent {width}x{int(width / max_aspect_ratio)} {preview_filename}",
                ignore_codes=[1],
            )
    
def preview(result: str, result_filename: str, schema: Union[str, None] = None) -> None:
    with open("result.tex", "w") as result_file:
        result_file.write(generate_tex(result, result_filename, schema)[0])
    render_preview("result.tex")

def get_webpage(address: str) -> str:
    response = urlopen(Request(address, headers={"User-Agent": "Mozilla/5.0"}))
    return response.read().decode()

def preview_from_url(url: str, schema: Union[str, None] = None) -> None:
    result_filename = url.split("/")[-1]
    result = get_webpage(url)
    preview(result_filename, result, schema)

'''
from IPython.display import Image

preview(
    "RNM_general__a_2_b_-1_-6_53642", 
    '[[[2, 11, 34, 46, 29, 7], [-1, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0], "6.943107048732200430941943037200926006848", [], [], 40]]'
)

Image("preview.png")
from IPython.display import Image

preview_from_url(
    "https://rnma.xyz/boinc/temp-data/output/RNM_PL5__a_2_b_-1_-6_53642",
    "general"
)

Image("preview.png")
'''