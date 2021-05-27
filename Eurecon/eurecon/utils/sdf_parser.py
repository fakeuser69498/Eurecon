import numpy as np

DEFAULT_SPAN = 10


class SdfDataObject:
    row_before_data = []
    row_after_data = []
    data_on_the_right = []
    data_on_the_left = []
    data = []
    type = None

    @classmethod
    def make_from_file(cls, input_directory):
        obj = cls()
        obj.data = []
        file = open(input_directory, "r")
        counter = 0
        ATOM_STARTED = False
        ATOMS_ENDED = False
        for row in file:
            counter += 1
            if "V2000" in row:
                obj.type = "V2000"

            if "V3000" in row:
                obj.type = "V3000"

            if counter > 4 and "." in row and not ATOMS_ENDED:
                ATOM_STARTED = True
                if obj.type == "V2000":
                    x = float(row[:10])
                    y = float(row[10:20])
                    z = float(row[20:30])
                    obj.data.append([x, y, z])
                    obj.data_on_the_right.append(row[30:])
                else:
                    x_span = obj.get_span()
                    x, y, z = row.split(" ")[5:8]
                    extra = row.split(" ")[8:]
                    obj.data.append([float(x), float(y), float(z)])
                    obj.data_on_the_left.append(row[:x_span])
                    obj.data_on_the_right.append(
                        obj.convert_data_from_list_to_row(extra)
                    )

            elif ATOMS_ENDED:
                obj.row_after_data.append(row)
                ATOMS_ENDED = True
            elif not ATOM_STARTED:
                obj.row_before_data.append(row)
                continue

        obj.data = np.asarray(obj.data)
        return obj

    def convert_data_from_list_to_row(self, list):
        result = ""
        for el in list:
            result = result + el + " "
        return result

    def get_span(self):
        current_mol_len = len(self.data) + 1
        margin = len(str(current_mol_len)) - 1

        return DEFAULT_SPAN + margin

    def configure_lines(self, data):
        result = []
        if self.type == "V2000":
            format = "{:>10.4f}{:>10.4f}{:>10.4f}"
            counter = 0
            for row in data:
                tmp = format.format(row[0], row[1], row[2])
                result.append(tmp + self.data_on_the_right[counter])
                counter += 1
        else:
            format_negative = "{:>8.4f}"
            format_positive = "{:>7.4f}"
            counter = 0
            for row in data:
                format = ""
                for el in row:
                    if 10.0 > el > 0.0:
                        format = format + format_positive
                    else:
                        format = format + format_negative
                tmp = format.format(row[0], row[1], row[2])
                result.append(
                    self.data_on_the_left[counter]
                    + tmp
                    + " "
                    + self.data_on_the_right[counter][:-1]
                )
                counter += 1

        return result

    def to_file(self, path):
        file = open(path, "w+")
        file.writelines(self.row_before_data)
        file.writelines(self.configure_lines(self.data))
        file.writelines(self.row_after_data)
