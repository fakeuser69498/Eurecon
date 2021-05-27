import os

try:
    from openbabel import pybel
except:
    raise Exception('OpenBabel import error')


def convert_input_file(path: str, format: str):
    """
        The function which converts the bio format file given as an input to .xyz extension.

        Parameters:
            path (String): Path to the bio format file.
            format (String): Given bio file extension.

        Returns:
            None.
    """
    file = next(pybel.readfile(format, path))
    file.write("xyz", path.replace(f".{format}", ".xyz"), True)


def convert_result_files(output_directory: str, to_format: str):
    """
        The function which converts the resulting conformations in .xyz format
            to an initial file extension.

        Parameters:
            output_directory (String): Path to the output directory.
            to_format (String): Resulting file extension.

        Returns:
            None.
    """
    result_files = os.listdir(output_directory)
    for file in result_files:
        if not file.endswith('.xyz'):
            continue
        xyz_file = next(pybel.readfile("xyz", f"{output_directory}/{file}"))
        result_file_name = file.replace(".xyz", f".{to_format}")
        xyz_file.write(to_format, f"{output_directory}/{result_file_name}", True)
        os.remove(f"{output_directory}/{file}")
