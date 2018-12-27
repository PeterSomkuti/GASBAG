import pandas as pd
import numpy as np


# Read in ephimeris data
df = pd.read_fwf('sun_ephem.txt')

# Make new column for day of year, but first, we have
# to convert the "MMM DD" strings into date time objects
date_list = []
for i in range(len(df)):
    date_list.append("2006-{:s}-{:d}".format(df['Mon'][i],
                                             df['Day'][i]))

df.index = pd.to_datetime(date_list)
df['DayOfYear'] = df.index.dayofyear

# Now we can perform a nice polynomial fit for both distance and
# velocity.

dist_fit = np.polyfit(df['DayOfYear'], df['Dist'], 8)
print("Distance polynomial coefficients [AU]: ")
for x in dist_fit[::-1]:
    print(x)

print("\n\n\n")

rv_fit = np.polyfit(df['DayOfYear'], df['RV'], 8)
print("Relative velocity polynomial coefficients [km/s]: ")
for x in rv_fit[::-1]:
    print(x)
