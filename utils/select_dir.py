import ROOT


def copy_directory(src_dir, dest_dir):
    # Loop over all keys in the source directory
    for key in src_dir.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            # If it's a directory, create a new one in the destination and recurse
            new_dir = dest_dir.mkdir(obj.GetName())
            copy_directory(obj, new_dir)
        else:
            # Otherwise, write the object to the destination directory
            dest_dir.cd()
            obj.Write()


def copy_tdirectory(input_file, src_dir_path, output_file, new_dir_name):
    # Open the input ROOT file
    infile = ROOT.TFile.Open(input_file, "READ")
    if not infile or infile.IsZombie():
        print(f"Cannot open input file: {input_file}")
        return

    # Get the source directory
    src_dir = infile.Get(src_dir_path)
    if not src_dir:
        print(f"Cannot find directory: {src_dir_path}")
        infile.Close()
        return

    # Open the output ROOT file
    outfile = ROOT.TFile.Open(output_file, "RECREATE")
    if not outfile or outfile.IsZombie():
        print(f"Cannot open output file: {output_file}")
        infile.Close()
        return

    # Create the new directory in the output file
    new_dir = outfile.mkdir(new_dir_name)
    copy_directory(src_dir, new_dir)

    # Write and close files
    outfile.Write()
    outfile.Close()
    infile.Close()
    print(f"Copied '{src_dir_path}' to '{output_file}:{new_dir_name}'")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 5:
        print(
            "Usage: python copy_tdirectory.py <input.root> <src_dir_path> <output.root> <new_dir_name>"
        )
        sys.exit(1)
    copy_tdirectory(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
