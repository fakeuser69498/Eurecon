import os.path as path
import sys
import os

import click


base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)

from eurecon import Eurecon




@click.command()
@click.option('--rmsd', '-r', help='Desired RMSD parameter')
@click.option('--partition', '-p', help='Desired partition parameter; 0 - translation only, 1 - rotation only')
@click.option('--axes_path', '-a', help='Path to the file, containing tessellation axes')
@click.option('--input_directory', '-i', help='Input directory for the base file')
@click.option('--output_directory', '-o', help='Output directory for the generated conformations')
@click.option('--stdout_mode', '-s', help='Enables/disables resulting file generation')
@click.option('--weights_file', '-w', default=None, help='Path to the file, containing weights')
@click.option('--debug', '-d', help='Enables resulting RMSD check for the generated conformations')
def main(rmsd, partition, axes_path, input_directory, output_directory, stdout_mode, weights_file, debug):
    """Hello."""
    path_list = os.listdir(input_directory)  # вместо refined-set указание папки (везде)
    result = []

    for path in path_list:
        try:
            for path2 in os.listdir(f'test_data/objects/{path}'):
                result.append(f'test_data/objects/{path}/{path2}')
        except Exception:
            pass

    data_for_plot = {}
    for object_path in result:
        print(object_path)
        format = object_path.split('.')[-1]
        eurecon = Eurecon(
            object_path,
            axes_path,
            output_directory,
            float(rmsd),
            float(partition),
            bool(int(stdout_mode)),
            weights_file,
            debug
        )

        _,_,_,metrics = eurecon.start()
        if format not in data_for_plot:
            data_for_plot[format] = {}

        data_for_plot[format][metrics['object_length']] = metrics['total_no_parsing']

    f = open('plot_one_random-klkllklklk.txt', 'w+')
    import json
    f.write(json.dumps(data_for_plot))


if __name__ == '__main__':
    main()
