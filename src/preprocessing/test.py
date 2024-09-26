import pandas as pd
import numpy as np
 
# data
Student = {
    'Name': ['John', 'Jay', 'sachin', 'Geetha', 'Amutha', 'ganesh'],
    'gender': ['male', 'male', 'male', 'female', 'female', 'male'],
    'math score': [50, 100, 70, 80, 75, 40],
    'test preparation': ['none', 'completed', 'none', 'completed',
                         'completed', 'none'],
}
 
# creating a Dataframe object
df = pd.DataFrame(Student)
 
print(df)
# Applying the condition
print(df.loc[df["Name"] == "Jay", "math score"].tolist()[0])