import os
import time

import inquirer
import pubchempy as pcp

from nomodeco.libraries import molpro_parser
from nomodeco.libraries import logfile


def _read_and_enumerate_xyz(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    atoms = []
    for i, line in enumerate(lines):
        if i >= 2:
            parts = line.split()
            if len(parts) == 4:
                atom, x, y, z = parts
                atoms.append((atom, float(x), float(y), float(z)))

    enumerated_atoms = [
        (f"{atom}{i}", x, y, z) for i, (atom, x, y, z) in enumerate(atoms, start=1)
    ]

    with open(file_path, "w") as f:
        f.write(f"{len(enumerated_atoms)} \n")
        f.write("File enumerated by Nomodeco.py \n")
        for atom, x, y, z in enumerated_atoms:
            f.write(f"{atom} {x: .4f} {y: .4f} {z: .4f}\n")


def _get_user_selection(choices):
    def _prompt():
        questions = [
            inquirer.Checkbox(
                "selected_files",
                message="Select structure for calculation (Spacebar to select)",
                choices=choices,
            ),
        ]
        answers = inquirer.prompt(questions)
        if not answers["selected_files"]:
            os.system("clear")
            print("Make at least one selection!")
            return _prompt()
        return answers["selected_files"]

    return _prompt()


def _build_isotope_string(change_to_isotope_lst):
    return "; ".join([f"mass, {v}=2.014" for v in change_to_isotope_lst]) + ";"


def _parse_output(output_file_path, change_to_isotope_lst=None):
    with open(output_file_path) as f:
        atoms = molpro_parser.parse_xyz_from_inputfile(f)
        n_atoms = len(atoms)
    if change_to_isotope_lst:
        for atom in atoms:
            if atom.symbol in change_to_isotope_lst:
                atom.swap_deuterium()
    with open(output_file_path) as f:
        CartesianF_Matrix = molpro_parser.parse_Cartesian_F_Matrix_from_inputfile(f)
        outputfile = logfile.create_filename_out(f.name)
    return atoms, n_atoms, CartesianF_Matrix, outputfile


def _run_s22():
    from ase.collections import s22
    from pymolpro import Project
    import ase.io

    time.sleep(1)
    while True:
        use_isotopes = input("Do you want to use isotopes in the calculation? [y/n] ")

        if use_isotopes == "n":
            print("Available Molecules:")
            structure_name = _get_user_selection(list(s22.names))[0]
            initial = s22[structure_name]

            p = Project(structure_name)
            ase.io.write(p.filename() + "/initial.xyz", initial)
            print("Now running molpro calculation ...")
            p.write_input(
                """
               orient,mass
               geometry=initial.xyz
               mass, iso
               basis=6-311g(d,p)
               {hf
               start, atden}
               optg;
               {frequencies, symm=auto, print=0, analytical}
               put, molden, %s.molden"""
                % structure_name
            )
            p.run(wait=True)
            assert p.status == "completed"
            print("... molpro calculation finished ...")
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path)

        if use_isotopes == "y":
            print("Available Molecules:")
            structure_name = _get_user_selection(list(s22.names))[0]
            initial = s22[structure_name]
            project_name = input("Specify Output Name: ")

            p = Project(project_name)
            ase.io.write(p.filename() + "/initial.xyz", initial)
            _read_and_enumerate_xyz(os.path.abspath(p.filename() + "/initial.xyz"))

            with open(os.path.abspath(p.filename() + "/initial.xyz")) as f:
                lines = [line.rstrip() for line in f]
            elements = [lines[i].strip()[:2] for i in range(2, len(lines))]
            print("The following atoms where found in the inputfile\n", elements)

            change_to_isotope_lst = input("Specify Atoms (use , as a seperator) ").split(",")
            isotope_string = _build_isotope_string(change_to_isotope_lst)

            print("Now running molpro calculation ...")
            p.write_input(
                """
                orient,mass
                geometry=initial.xyz
                %s ! Set custom mass for hydrogen
                basis=6-311g(d,p)
                {hf
                start, atden}
                optg;
                {frequencies, symm=auto, print=0, analytical}
                put, molden, %s.molden"""
                % (isotope_string, project_name)
            )
            p.run(wait=True)
            assert p.status == "completed"
            print("... molpro calculation finished ...")
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path, change_to_isotope_lst)


def _run_g2():
    from ase.collections import g2
    from pymolpro import Project
    import ase.io

    time.sleep(1)
    while True:
        use_isotopes = input("Do you want to use isotopes in the calculation? [y/n] ")

        if use_isotopes == "n":
            print("Available Molecules:")
            structure_name = _get_user_selection(list(g2.names))[0]
            initial = g2[structure_name]

            p = Project(structure_name)
            ase.io.write(p.filename() + "/initial.xyz", initial)
            print("Now running molpro calculation ...")
            p.write_input(
                """
                orient,mass
                geometry=initial.xyz
                mass, iso
                basis=6-311g(d,p)
                {hf
                start, atden}
                optg;
                {frequencies, symm=auto, print=0, analytical}
                put, molden, %s.molden"""
                % structure_name
            )
            p.run(wait=True)
            assert p.status == "completed"
            print("... molpro calculation finished ...")
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path)

        if use_isotopes == "y":
            print("Available Molecules:")
            structure_name = _get_user_selection(list(g2.names))[0]
            initial = g2[structure_name]
            project_name = input("Specify Output Name: ")

            p = Project(project_name)
            ase.io.write(p.filename() + "/initial.xyz", initial)
            _read_and_enumerate_xyz(os.path.abspath(p.filename() + "/initial.xyz"))

            with open(os.path.abspath(p.filename() + "/initial.xyz")) as f:
                lines = [line.rstrip() for line in f]
            elements = [lines[i].strip()[:2] for i in range(2, len(lines))]
            print("The following atoms where found in the inputfile\n", elements)

            change_to_isotope_lst = input("Specify Atoms (use , as a seperator) ").split(",")
            isotope_string = _build_isotope_string(change_to_isotope_lst)

            print("Now running molpro calculation ...")
            p.write_input(
                """
                orient,mass
                geometry=initial.xyz
                %s ! Set custom mass for hydrogen
                basis=6-311g(d,p)
                {hf
                start, atden}
                optg;
                {frequencies, symm=auto, print=0, analytical}
                put, molden, %s.molden"""
                % (isotope_string, project_name)
            )
            p.run(wait=True)
            assert p.status == "completed"
            print("... molpro calculation finished ...")
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path, change_to_isotope_lst)


def _run_xyz_file():
    from pymolpro import Project

    while True:
        use_isotopes = input("Do you want to use isotopes in the calculation? [y/n] ")

        if use_isotopes == "n":
            structure_name = input("Specify Output Name: ")
            print("Following Listing files in the current directory ...")
            time.sleep(1)
            files = [f for f in os.listdir(".") if os.path.isfile(f) and f.endswith(".xyz")]
            questions = [
                inquirer.List(
                    ".xyz",
                    message="Following .xyz files are available in the working directory",
                    choices=files,
                ),
            ]
            xyz_file = inquirer.prompt(questions)[".xyz"]
            xyz_abs_path = os.path.abspath(xyz_file)

            print("Now running molpro calculation ...")
            p = Project(structure_name)
            p.write_input(
                """
              orient,mass
              geometry=%s
              mass, iso
              basis=6-311g(d,p)
              {hf
              start, atden}
              optg;
              {frequencies, symm=auto, print=0, analytical}
              put, molden, %s.molden"""
                % (xyz_abs_path, structure_name)
            )
            p.run(wait=True)
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path)

        elif use_isotopes == "y":
            structure_name = input("Specify Output Name: ")
            print("Following Listing files in the current directory ...")
            time.sleep(1)
            files = [f for f in os.listdir(".") if os.path.isfile(f) and f.endswith(".xyz")]
            questions = [
                inquirer.List(
                    ".xyz",
                    message="Following .xyz files are available in the working directory",
                    choices=files,
                ),
            ]
            xyz_file = inquirer.prompt(questions)[".xyz"]

            with open(xyz_file) as f:
                lines = [line.rstrip() for line in f]
            elements = [lines[i].strip()[:2] for i in range(2, len(lines))]
            print("The following atoms where found in the inputfile\n", elements)

            change_to_isotope_lst = input("Specify Atoms (use , as a seperator) ").split(",")
            isotope_string = _build_isotope_string(change_to_isotope_lst)
            xyz_abs_path = os.path.abspath(xyz_file)

            print("Now running molpro calculation ...")
            p = Project(structure_name)
            p.write_input(
                """
             orient,mass
             geometry=%s
             %s ! Set custom mass for hydrogen
             basis=6-311g(d,p)
             {hf
             start, atden}
             optg;
             {frequencies, symm=auto, print=0, analytical}
             put, molden, %s.molden"""
                % (xyz_abs_path, isotope_string, structure_name)
            )
            p.run(wait=True)
            os.environ["OUT_FILE_LINK"] = p.output_file_path
            return _parse_output(p.output_file_path, change_to_isotope_lst)


def _run_pubchem():
    from pymolpro import Project

    os.system("clear")
    user_search_input = input("Specify molecule to search for: ")
    compounds = pcp.get_compounds(user_search_input, "name")

    molecule_form_iso_smiles = []
    for i, compound in enumerate(compounds):
        molecule_form_iso_smiles.append({
            "Idx": i,
            "Chemical_Formula": compound.molecular_formula,
            "Isomeric_Smiles": compound.isomeric_smiles,
            "Pubchem_CID": compound.cid,
        })

    def get_user_selection():
        questions = [
            inquirer.Checkbox(
                "selected_files",
                message="Following Compounds where found (Spacebar to select):",
                choices=molecule_form_iso_smiles,
            ),
        ]
        answers = inquirer.prompt(questions)
        if not answers["selected_files"]:
            os.system("clear")
            print("Make at least one selection!")
            return get_user_selection()
        return answers["selected_files"]

    compound_cid = get_user_selection()[0]["Pubchem_CID"]
    molecule_selected = pcp.Compound.from_cid(compound_cid, record_type="3d")
    dict_atom_coords = molecule_selected.to_dict(properties=["atoms"])

    xyz_filename = user_search_input + ".xyz"
    with open(xyz_filename, "w") as f:
        f.write(f"{len(dict_atom_coords['atoms'])}\n")
        f.write("XYZ File generated by Nomodeco.py\n")
        for atom in dict_atom_coords["atoms"]:
            f.write(f"{atom['element']}{atom['aid']} {atom['x']:.4f} {atom['y']:.4f} {atom['z']:.4f}\n")

    xyz_abs_path = os.path.abspath(xyz_filename)
    print("Now running molpro calculation ...")
    p = Project(user_search_input)
    p.write_input(
        """
        orient,mass
        geometry=%s
        mass, iso
        basis=6-311g(d,p)
        {hf
        start, atden}
        optg;
        {frequencies, symm=auto, print=0, analytical}
        put, molden, %s.molden"""
        % (xyz_abs_path, user_search_input)
    )
    p.run(wait=True)
    os.environ["OUT_FILE_LINK"] = p.output_file_path
    return _parse_output(p.output_file_path)


def run_pymolpro_workflow():
    """Run the interactive pymolpro calculation workflow.

    Returns:
        tuple: (atoms, n_atoms, CartesianF_Matrix, outputfile)
    """
    questions = [
        inquirer.List(
            "mode",
            message="Select calculation mode",
            choices=[
                "Use g2 databank",
                "Use s22 databank",
                "Do pubchem search",
                ".xyz file calculation",
            ],
        ),
    ]
    answers = inquirer.prompt(questions)
    mode = answers["mode"]

    if mode == "Use s22 databank":
        return _run_s22()
    elif mode == "Use g2 databank":
        return _run_g2()
    elif mode == ".xyz file calculation":
        return _run_xyz_file()
    elif mode == "Do pubchem search":
        return _run_pubchem()