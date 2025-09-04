import os, sys, argparse

def repair_outdock(file_path, output_path):
    with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            tokens = line.strip().split()

            if len(tokens) == 20 and 'ZINC' in tokens[1]:
                # Calculate score
                terms = [float(x) for x in tokens[11:16]]
                score = round(sum(terms), 2)

                # Clean up the last element
                fixed_rhyd = tokens[-1].replace('*', '')

                # Find where last element begins in the original line
                last_start = line.rfind(tokens[-1])

                # Prefix = everything before last element (with original spacing preserved)
                prefix = line[:last_start]

                # Build new line: prefix + cleaned last col + tab + score
                new_line = f"{prefix}{fixed_rhyd}\t{score}\n"

                outfile.write(new_line)
            else:
                outfile.write(line)

def main():

    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('file_input_path') 
    parser.add_argument('file_save_path') 
    args = parser.parse_args()

    repair_outdock(args.file_input_path, args.file_save_path)


if __name__ == '__main__':
    main()
    