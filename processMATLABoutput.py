import csv

def convert_txt_to_csv(open_path, save_path, excluded_variables):
    data = {}
    with open(open_path, 'r') as file:
        # Skip the first three lines
        for _ in range(3):
            next(file)

        for line in file:
            # Assuming each line is in the format "Variable | Value"
            parts = line.split('|')
            if len(parts) == 2:
                variable = parts[0].strip()
                value = parts[1].strip()
                if variable not in excluded_variables:
                    data[variable] = value

    with open(save_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Variable', 'Value'])
        for variable, value in data.items():
            writer.writerow([variable, value])

# Example usage
open_path = '/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs/mainvariables20240113214128.txt'  # Replace with your text file path
save_path = '/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs/matlabmainvariables.csv'  # Replace with your desired CSV file path
excluded_variables = ['filename']  # Replace with variables to exclude

convert_txt_to_csv(open_path, save_path, excluded_variables)
