[![DOI](https://zenodo.org/badge/971190779.svg)](https://doi.org/10.5281/zenodo.15309368) [![Launch App](https://img.shields.io/badge/Launch%20App-Shiny-blue?logo=R)](https://anatomeshr.serve.scilifelab.se/)

# Anatomeshr

> _A converter of 2D plant anatomical data into finite element mesh (GEO) format_

**Lead development:**  Adrien Heymans

**Contributors:** Ioannis Theodorou, Grégoire Loupit, Vinod Kumar, Gonzalo Revilla, Olivier Ali, Stéphanie Robert, Stéphane Verger

---

## 1. About

[**Anatomeshr**](https://anatomeshr.serve.scilifelab.se/) is a user-friendly **Shiny app** that converts 2D plant tissue image data into **finite element-compatible meshes** for biomechanical simulations.

This tool bridges high-resolution anatomical imaging (CellSeT, ImageJ) and simulation frameworks like [**BVPy**](https://mosaic.gitlabpages.inria.fr/bvpy/), providing:

- Input support for:
  - `CellSeT` or [`GRANAR`](https://granar.github.io/) `.xml` files
  - `ImageJ` `.roi` files (single cells)
  - Predefined anatomical templates
- Output: `.geo` files with:
  - Customizable wall thickness
  - Adjustable vertex smoothing
  - Cell Tags to easily mark boundary conditions 
- Seamless integration with finite element software (e.g., Fenics/BVPy)

[**Start the Shiny app**](https://anatomeshr.serve.scilifelab.se/)

## 2. Installation

### Clone the repository

```bash
git clone https://github.com/SRobertGroup/Anatomeshr.git
cd Anatomeshr
```

### R Dependencies

Ensure you're using **R 4.3+**, then install required packages:

```r
install.packages(c("sf", "units", "xml2", "tidyverse", "viridis",
                   "smoothr", "RImageJROI", "utils", "codetools", "purrr"),
                 dependencies = TRUE)
```

### Run via Docker

To run with Docker (recommended):

```bash
docker pull heymansadrien/anatomeshr:0.0.6
```

Then launch it:

```bash
docker run -p 3838:3838 heymansadrien/anatomeshr:0.0.6
```

App will be available at: [http://localhost:3838](http://localhost:3838)

## 3. Usage

### Step 1: Upload Data

Choose one of the input types:

- 📂 **CellSet XML:** `.xml` exported from [CellSet](https://www.cellset.org/) or GRANAR  
- 📂 **ImageJ ROIs:** One or more `.roi` files (drawn manually in ImageJ/Fiji)
- 📂 **preloaded anatomies** from literature

### Step 2: Prepare Geometry

Customize geometry generation:

- Adjust **cell wall thickness**
- Tweak **corner smoothing** for more realistic cell shapes

### Step 3: Export Mesh

Export your mesh as:

- `plot.png`, `plot.svg` — anatomical visualization
- `geometry.csv` — tabular coordinate of the geo data
- `geometry.geo` — mesh input file for finite element analysis through [GMSH](https://gmsh.info/) 


> ⚠️ **Tips & Warnings**
>
> - Concave cells in `.xml` input can cause meshing errors
> - When using `.roi` input, consider increasing wall thickness to improve polygon closure
> - Use higher smoothing values to improve mesh quality in large or irregular cells


## 4. Citation

> Please cite this repository if you use **Anatomeshr** in your work.  
> A manuscript describing the method is currently in preparation.

---

## Acknowledgments

Anatomeshr was developed by members of the [**Stéphanie Robert Group**](https://srobertgroup.com/) and [**Stéphane Verger Group**](https://www.upsc.se/researchers/6177-verger-stephane-mechanics-and-dynamics-of-cell-to-cell-adhesion-in-plants.html), in collaboration with plant mechanobiology researchers.
