from pylatex import Document, Section, Alignat


# Create a LaTeX document with a list of equations
def generate_latex(file_name, eqns=None):
    if eqns is None:
        eqns = []
    doc = Document()

    with doc.create(Section("Automatic Conjectures")):
        doc.append("These are the conjectures detected by the algorithm.")

        for eqn in eqns:
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(eqn)

    doc.generate_tex(file_name)
    # doc.generate_pdf(file_name, clean_tex=False)
