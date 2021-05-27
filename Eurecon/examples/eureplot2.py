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
    tess_path_list = os.listdir(axes_path)
    tess_result = []

    for path in tess_path_list:
        try:
            tess_result.append(f'test_data/tess/{path}')
        except Exception:
            pass

    path_list = os.listdir(input_directory)
    obj_result = []

    for path in path_list:
        try:
            for path2 in os.listdir(f'test_data/objects/{path}'):
                obj_result.append(f'test_data/objects/{path}/{path2}')
        except Exception:
            pass

    some_res = []
    for el in obj_result:
        if '96.' in el:
            some_res.append(el)

    data_for_plot = {}
    for object_path in some_res:
        print(object_path)
        for axes_path in tess_result:
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

            if object_path.endswith('_.ply'):
                print('!')
                continue

            if format not in data_for_plot:
                data_for_plot[format] = {}

            data_for_plot[format][metrics['axes_amount']] = metrics['total_no_parsing']

    f = open('plot_two.txt', 'w+')
    import json

    f.write(json.dumps(data_for_plot))
    f.close()
if __name__ == '__main__':
    main()
