# Midnight-Club-3-PKG-Visualizer by JD3 Ent and MC3x

## Overview
This tool allows users to view and export 3D models from Midnight Club 3 `.pck` files. 
It parses vertex, UV, and face data, organizes them into structured meshes, and displays them in an interactive 3D viewer using PyVista and PyQt5.
It will receive updates and changes but this is the base code that can be used by MC3x, and I (JD3 Ent) will be contributing to this project on my free time.
My first commit will be english translation.
---

## Features
- Parse PCK files to extract 3D model data  
- Detect and group vertices into proper submeshes  
- Identify active and inactive faces for proper rendering  
- Display UV mapping information  
- Interactive 3D visualization with zoom, rotate, and pan controls  
- Customization settings for wireframe mode, shading, transparency, and background color  
- Toggle visibility of inactive faces  
- Export models to OBJ format for use in other 3D software  

---

## How It Works

### Step 1: Open a PCK File
1. Launch the program.
2. Click **"Abrir PCK"** to select a `.pck` file.
3. The file is parsed, and vertices, UVs, and faces are extracted.

### Step 2: View the Model
- The extracted vertex groups are displayed in different colors.
- Active faces are rendered normally, while inactive faces can be toggled on or off.
- The Information Panel provides details about each mesh.

### Step 3: Customize the View
- Adjust settings via the **Settings Panel**:  
  - Change background color  
  - Toggle wireframe mode  
  - Enable or disable grid and axes  
  - Adjust shading (Flat or Smooth)  
  - Modify transparency  
  - Scale the model by dividing vertex coordinates by 256  

### Step 4: Export the Model
- Click **"Exportar para OBJ"** to save the model in `.obj` format.
- UV mapping is preserved in the exported file.

---

## Technical Details

### File Parsing (`PCKParser`)
The parser detects vertex groups using specific binary patterns:

- **Vertex headers:**  
  - `EE 00 XX 69`  
  - `1B 02 XX 69`  
- **UV headers:**  
  - `C4 00 XX 65`  
  - `F1 01 XX 65`  
- **Face data headers:**  
  - `9A 00 XX 6A`  
  - `C7 01 XX 6A`  

### Data Extraction Process
1. Read vertex data (each vertex is stored as `(x, y, z)` 16-bit integers).  
2. Extract UV coordinates (each UV pair is stored as `(u, v)` 16-bit signed values).  
3. Identify faces and classify them as active or inactive.  
4. Convert local indices to global indices for correct mesh assembly.  

### 3D Rendering (`MainWindow`)
- Uses PyVista to visualize the extracted model.  
- Supports real-time modifications (wireframe mode, shading, transparency).  
- Renders each vertex group as a separate mesh with assigned colors.  
- Provides an interactive UI with model settings and export options.  

---

## Controls & UI Panels

| Panel                | Functionality |
|----------------------|--------------------------------|
| **Information Panel** | Displays vertex, UV, and face counts for each mesh. |
| **Toggle Panel**     | Allows users to toggle inactive face visibility. |
| **Settings Panel**   | Controls shading, background color, wireframe mode, and transparency. |

---

## Installation & Dependencies

### Requirements
Ensure you have Python 3 installed and install the required dependencies:
