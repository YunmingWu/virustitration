# Virus Titer Analysis App

This is a Shiny app designed for the analysis of quantitative polymerase chain reaction (qPCR) data. 
It helps researchers perform a linear regression on their standard curve data, calculate key metrics like PCR efficiency, and predict the virus titer of unknown samples.
This app is based on the protocol describe here:


## How to use

### Prerequisites:
You need to have R and RStudio installed on your computer.
The app also requires the following R packages: shiny ggplot2 DT.
The app has a built-in check to automatically install any missing packages. However, you can manually install them by running this command in your R console:
install.packages(c("shiny", "ggplot2", "DT"))

### Download:
First download the entire folder to your local computer.

### Run the App: 
Open the app.R file in RStudio and click the "Run App" button at the top-right of the script editor.

### In the app:
#### Step I: Input Your Standard Curve Data
An editable table where you can enter the Cq and SQ values from your standard curve.

#### Step II: Predict Sample Quantity (SQ)
Enter the Cq value of your unknown sample.
The app will display the predicted starting quantity (SQ) of that sample.

#### Step III: Calculate Virus Titer
Enter the viral genome length (in base pairs) and the dilution factor.
The app will display the final calculated virus titer in "genomic copies/mL".



## How It Works
The app uses the following formulas for its calculations:

Standard Curve Equation:
y=mx+b
where y is the Cq value and x is the log$_{10}$(SQ).


Virus Titer:
The formula used for this calculation is:
Titer = {Sq * 10-9 * 10^23 (molar to molecule number) * dilution factor} / { genome plasmid length * 650 (average molecular weight of a single nucleotide) }

SQ: The predicted starting quantity of the sample.


