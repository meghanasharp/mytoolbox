# Define a function that will read a generic text file and pull any keys that have the format:
#   Keyword: Value
#     NOTE: The ": " colon/space sequence is the only criterion needed to specify the name.
#             E.G. Input_file: C:\myfile.shp
#           Lines without the Keyword: Value sequence will be ignored.
#   Returns: Dictionary of key/value pairs

def ReadParamFileGeneric(paramfile):
    # Test if the file exists.
    if not os.path.isfile(paramfile):
        print("This file not found: " + paramfile)
        return False, None

    # Set up a blank dictionary value.
    params = {}

    # Try opening the file, using the "with ... as ..." construct, which
    #   is a wrapper around the file I/O to ensure that the file
    #   tunnel always gets cleaned up, even when there's an error.
    try:
        with open(paramfile) as f:
            # Go through the file line by line and split each line on the colon/space combo
            for line in f:
                check = line.split(": ")
                # If the result has exactly two values, we know we've found one match. Add this to the dictionary.
                if len(check) == 2:
                    # By assigning a new key/value pair, the dictionary grows.
                    # We also need to make sure that any leading/extra spaces in the value are cleaned up with .strip().
                    params[check[0]] = check[1].strip()
                    
                    # If the length is not exactly 2, we assume this is not an interesting line.
                    #   This could cause problems if you want to pass multiple things on a single line.
                    #   How could you adapt this to be more flexible?
                    #     E.G. elif len(check) > 2: ... for i, kv in enumerate(check): ... etc.
                    
    except:  # We get here if it couldn't open the file
        print("Problem reading the paramfile: " + paramfile)
        return False, None

    # If we got here (note the indentation), we made it through the loop and can return the params dictionary object.
    # Return the dictionary value, along with a success flag.
    return True, params