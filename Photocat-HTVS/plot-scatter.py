import sys,os
import pandas as pd
import matplotlib.pyplot as plt

title='%s' % sys.argv[1].split('.')[0]

# Step 1: Read the data from the CSV file
csv_name=title+'.csv'
df = pd.read_csv(csv_name)

# Step 2: Extract the data for the scatter plot
xvalue='Oxi_'+title
yvalue='Red_'+title
x_values = df[xvalue]
y_values = df[yvalue]

# Step 3: Change the font settings
font = {'family': 'sans-serif', 'weight': 'normal', 'size': 14}
plt.rc('font', **font)
#plt.figure(figsize=(8, 4))

# Step 4: Create the scatter plot with colormap
color_hex='#154c79'
plt.scatter(x_values, y_values,c=color_hex)
plt.xlabel('Oxidation Potential (eV)')
plt.ylabel('Reduction Potential (eV)')
plt.grid(False)

# Step 6: Save the plot as a PNG image
plt.savefig('plot_'+title+'.png',transparent=True, bbox_inches='tight',dpi=850)

