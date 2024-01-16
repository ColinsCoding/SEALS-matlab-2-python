import csv

def convert_txt_to_csv(open_path, save_path, excluded_variables):
    data = []
    with open(open_path, 'r') as file:
        # Skip the first three lines
        for _ in range(3):
            next(file)

        for line in file:
            # Now each line is in the format "Variable | Value | Datatype"
            parts = line.split('|')
            if len(parts) == 3:
                variable = parts[0].strip()
                value = parts[1].strip()
                datatype = parts[2].strip()
                if variable not in excluded_variables:
                    data.append([variable, value, datatype])

    with open(save_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Variable', 'Value', 'Datatype'])
        for row in data:
            writer.writerow(row)


# Example usage
open_path = '/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs/matlabmainvariables20240115115635.txt'  # Replace with your text file path
save_path = '/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs/matlabmainvariables.csv'  # Replace with your desired CSV file path
excluded_variables = ['filename']  # Replace with variables to exclude

convert_txt_to_csv(open_path, save_path, excluded_variables)
