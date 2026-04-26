# Project README

## Overview
This project is a simple 3D engine written in C. It allows for the creation and display of basic 3D scenes using a wireframe rendering technique.

## Features
- Basic 3D rendering
- Wireframe drawing
- Simple setup and execution

## Project Structure
```
<Gui_3D_Engine_2/>
├── build/
├── src/
│   ├── Main.c
│   └── *.h
├── Makefile.linux
├── Makefile.windows
├── Makefile.wine
└── README.md
```

### Prerequisites
- C/C++ Compiler and Debugger (GCC, Clang)
- Make utility
- Standard development tools
- X11 for Linux

## Build & Run
To build the project on Linux:

```bash
cd <Gui_3D_Engine_2>
make -f Makefile.linux all
```

To run the project:

```bash
make -f Makefile.linux exe
```

For Windows, ensure you have MinGW installed and then:

```bash
cd <Gui_3D_Engine_2>
make -f Makefile.windows all
make -f Makefile.windows exe
```

For Wine on Linux:

```bash
cd <Gui_3D_Engine_2>
make -f Makefile.wine all
make -f Makefile.wine exe
```

To build for web, ensure you have Emscripten installed and then:

```bash
cd <Gui_3D_Engine_2>
make -f Makefile.web all
make -f Makefile.web exe
```

Note: The project assumes a basic understanding of C programming and the use of make utilities for building projects.