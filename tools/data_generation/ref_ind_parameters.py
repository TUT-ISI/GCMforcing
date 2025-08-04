import pandas as pd

"""
by  Atte Laakso / Aalto University

Reads refractive indices of chemical compounds for different models
"""


# Function for saving data from a csv file to a dictionary
def load_data(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path, delimiter=';', index_col=0)
    
    # The first row contains the compound names
    names_of_compounds = df.columns
    
    # The first column contains the model names
    models = df.index
    
    # Create a dictionary to store the data
    data = {}

    # Iterate over each row in the dataframe
    for index, row in df.iterrows():
        if index in ('EC-Earth3-AerChem'):
            print('Ref ind iteration skipped duo to lack of data')
            continue
        else:
            # index = model's name
            # Create a dictionary to store the data from this model
            modeldata = {}
            
            for compounds in names_of_compounds:
                # Intitialise a refrac dictionary for models refractive indices
                refrac = {}            
                # Extract real and imaginary parts
                value = row[compounds]
                # If not already, make the value a string
                if not isinstance(value, str):
                    value=str(value)
                # Reformate the value to mach the needs of the output
                value=value.lower().replace(' ', '').replace('i', '')
                if value == "":
                    continue
                # Check if splitting the value is possible
                parts=value.strip('()').split('+')
                if len(parts) != 2:
                    print('couldnt handle number'+str(value))
                    continue
                # Form the complex number (refractive index, absorption coefficient)
                real_part, imag_part = map(float, parts)
                modeldata[compounds] = complex(real_part, imag_part)
                    
            # Save the corresponding refractive indices under the name of the model      
            data[index] = modeldata
    return data

# Method for searching model data from the csv file containing refractive indices
def get_ref_ind(model_name, data):
    search_name=''
    # Search model name from the indices data
    for key in data.keys():
        if key.lower() in model_name.lower():
            search_name=key
    # Return data of given model if possible
    return data.get(search_name, "Model not found")

# Method to be called from elsewhere
def import_refrac(name):
    # Define file path for refractive indices storage
    file_path='../tools/data_generation/ref_ind_storage.csv'
    # Read refractive index data
    data=load_data(file_path)
    # Search refractive indices for given model and return those indices
    refrac=get_ref_ind(name, data)    
    return refrac



