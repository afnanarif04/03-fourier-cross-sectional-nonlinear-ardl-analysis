# Data Source

The panel dataset in `data.xlsx` was constructed from publicly available data.

## Primary source

**Our World in Data — CO2 and Greenhouse Gas Emissions / Energy**

- Dataset page: <https://github.com/owid/co2-data>
- Direct download (CSV): <https://nyc3.digitaloceanspaces.com/owid-public/data/co2/owid-co2-data.csv>
- Codebook (variable definitions): <https://github.com/owid/co2-data/blob/master/owid-co2-codebook.csv>

The Our World in Data energy and emissions data are compiled from the Energy Institute Statistical Review of World Energy, the U.S. Energy Information Administration, and the Global Carbon Project.

## Coverage used in this study

- Cross-sectional units: [e.g., 37 OECD countries]
- Time span: [e.g., 1990–2022]
- Variables retained: country identifier, year, and the two series used as `y` and `x`

## Construction notes

The raw download was filtered to the units and years above, missing observations were removed, and the series were arranged into the four-column panel format (`id`, `time`, `y`, `x`) used by the estimation script. No transformation other than the standard logarithm was applied to the level variables.

> Note: replace the bracketed placeholders above with the exact units, years, and variables used, and confirm the direct download link is current. Our World in Data occasionally updates the file path.
