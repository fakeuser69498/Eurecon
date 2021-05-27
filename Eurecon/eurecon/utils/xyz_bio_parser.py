import numpy as np

class XyzBioDataObject:
    row_before_data = []
    data_on_the_left = []
    data = []

    @classmethod
    def make_from_file(cls, input_directory):
        obj = cls()
        obj.data = []
        file = open(input_directory, "r")
        counter = 0
        for row in file:
            counter += 1
            if counter <= 2:
                obj.row_before_data.append(row)
                continue

            row = row.split()
            obj.data_on_the_left.append(row[0])

            x = float(row[1])
            y = float(row[2])
            z = float(row[3])
            obj.data.append([x, y, z])

        obj.data = np.asarray(obj.data).T
        return obj

    def configure_lines(self):
        result = []
        format = "{}{:>17.5f}{:>15.5f}{:>15.5f}\n"
        counter = 0
        for row in self.data.T:
            tmp = format.format(self.data_on_the_left[counter], row[0], row[1], row[2])
            result.append(tmp)
            counter += 1

        return result

    def to_file(self, path):
        file = open(path, "w+")
        file.writelines(self.row_before_data)
        file.writelines(self.configure_lines())


def is_bio(file_name: str) -> bool:
    bio_file = False
    with open(file_name) as file:
        for i, line in enumerate(file):
            if i == 3:
                data = line.split()
                try:
                    first_el = float(data[0])
                except ValueError:
                    bio_file = True
                break

    return bio_file
