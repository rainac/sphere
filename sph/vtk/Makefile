
sph_root = ..

test_base := test-cases
test_stem := test
test_file := $(test_stem).vtk
xml_suffix := .sphscene.xml

sph_base := ..
scripts  := $(sph_base)/sh/scenes

dxKlein2  := 0.01
dxMittel2 := 0.005
dxGross2  := 0.0025
dxRiesig2 := 0.00125

tests := \
1D/stempel/klein 1D/stempel/mittel 1D/stempel/gross \
2D/bar/klein/V001 2D/bar/klein/V010 2D/bar/klein/V100 \
2D/bar/mittel/V001 2D/bar/mittel/V010 2D/bar/mittel/V100 \
2D/bar/gross/V001 2D/bar/gross/V010 2D/bar/gross/V100 \
2D/line/mittel/V001 2D/line/mittel/V010 2D/line/mittel/V100 \
2D/monaghan-ellipsis/klein 2D/monaghan-ellipsis/mittel 2D/monaghan-ellipsis/gross \
3D/monaghan-lentil/klein 3D/monaghan-lentil/mittel 3D/monaghan-lentil/gross \
3D/monaghan-cigar/klein 3D/monaghan-cigar/mittel 3D/monaghan-cigar/gross \
2D/table-shot/klein/V001 2D/table-shot/klein/V010 2D/table-shot/klein/V100 \
2D/table-shot/mittel/V001 2D/table-shot/mittel/V010 2D/table-shot/mittel/V100 \
2D/table-shot/gross/V001 2D/table-shot/gross/V010 2D/table-shot/gross/V100 \
3D/table-shot/klein/V001 3D/table-shot/klein/V010 3D/table-shot/klein/V100 \
3D/table-shot/mittel/V001 3D/table-shot/mittel/V010 3D/table-shot/mittel/V100 \
2D/fountain/klein 2D/fountain/mittel 2D/fountain/gross \
3D/fountain/klein 3D/fountain/mittel 3D/fountain/gross \
2D/dambreak/klein 2D/dambreak/mittel 2D/dambreak/gross 2D/dambreak/riesig \
3D/dambreak/klein/216 3D/dambreak/klein/1000 3D/dambreak/klein/2744 \
3D/dambreak/mittel/10000 3D/dambreak/mittel/50000 \
3D/dambreak/gross/1000000 \
1D/burgers/klein \
2D/poiseuille-flow/klein 2D/poiseuille-flow/mittel 2D/poiseuille-flow/gross \
2D/poiseuille-flow-hex/klein 2D/poiseuille-flow-hex/mittel 2D/poiseuille-flow-hex/gross \
2D/falling-film/klein 2D/falling-film/mittel 2D/falling-film/gross \
2D/falling-film-vert/klein 2D/falling-film-vert/mittel 2D/falling-film-vert/gross \
2D/falling-film-triang/klein 2D/falling-film-triang/mittel 2D/falling-film-triang/gross \
2D/falling-film-hex/klein 2D/falling-film-hex/mittel 2D/falling-film-hex/gross \
2D/falling-film-2/klein 2D/falling-film-2/mittel 2D/falling-film-2/gross \
2D/falling-film-frei/klein 2D/falling-film-frei/mittel 2D/falling-film-frei/gross \
3D/falling-film/klein 3D/falling-film/mittel 3D/falling-film/gross \
2D/treppe3/klein 2D/treppe3/mittel 2D/treppe3/gross \
2D/treppe5/klein 2D/treppe5/mittel 2D/treppe5/gross \
2D/treppe11/klein 2D/treppe11/mittel 2D/treppe11/gross \
2D/treppe12/klein 2D/treppe12/mittel 2D/treppe12/gross \
2D/treppe21/klein 2D/treppe21/mittel 2D/treppe21/gross \
3D/treppe3/klein 3D/treppe3/mittel 3D/treppe3/gross \
3D/treppe5/klein 3D/treppe5/mittel 3D/treppe5/gross \
2D/box-drop/klein 2D/box-drop/mittel 2D/box-drop/gross \
2D/box/klein 2D/box/mittel 2D/box/gross \
3D/box-drop/klein 3D/box-drop/mittel 3D/box-drop/gross \
2D/utube/klein 2D/utube/mittel 2D/utube/gross \
2D/artesian-fountain/klein 2D/artesian-fountain/mittel 2D/artesian-fountain/gross 



#2D/box-filled/klein 2D/box-filled/mittel 2D/box-filled/gross \
#2D/box-filled-offset/klein 2D/box-filled-offset/mittel 2D/box-filled-offset/gross \
#2D/dambreak-obstacle/klein 2D/dambreak-obstacle/mittel 2D/dambreak-obstacle/gross# 2D/dambreak-obstacle/riesig \
#3D/dambreak-obstacle/klein 3D/dambreak-obstacle/mittel 3D/dambreak-obstacle/gross \

# 3D/dambreak/gross/200000 \
# 3D/dambreak/gross/2000000 \

#3D/pattern-test/klein 3D/pattern-test/mittel 3D/pattern-test/gross \
#2D/pattern-test/klein 2D/pattern-test/mittel 2D/pattern-test/gross \
#2D/pattern-rect-cart 2D/pattern-rect-hex 2D/pattern-rect-hex-odd 

# 3D/dambreak-obstacle/riesig

# 2D/table-column/klein/H10 2D/table-column/klein/H20 2D/table-column/klein/H40 \
# 2D/table-column/mittel/H10 2D/table-column/mittel/H20 2D/table-column/mittel/H40 \
# 3D/table-column/klein/H10 3D/table-column/klein/H20 3D/table-column/klein/H40 \
# 3D/table-column/mittel/H10 3D/table-column/mittel/H20 3D/table-column/mittel/H40 \
# klein/3D/table/125 klein/3D/table/27 mittel/3D/table/10000 mittel/3D/table/50000 \
# klein/3D/table-skewed/125 klein/3D/table-skewed/27 mittel/3D/table-skewed/10000 mittel/3D/table-skewed/50000 \
# mittel/2D/box-filled/25 mittel/2D/box-filled/121 mittel/2D/box-filled/625 mittel/2D/box-filled/2000 \
# mittel/2D/box-drop/25 mittel/2D/box-drop/121 mittel/2D/box-drop/625 mittel/2D/box-drop/2000 \
# mittel/2D/box-filled-large/121 mittel/2D/box-filled-large/625 mittel/2D/box-filled-large/2000 mittel/2D/box-filled-large/10000 \
# mittel/2D/box-drop-large/121 mittel/2D/box-drop-large/625 mittel/2D/box-drop-large/2000 mittel/2D/box-drop-large/10000 \
# mittel/3D/box-drop-large/261 mittel/3D/box-drop-large/1728 mittel/3D/box-drop-large/8000 mittel/3D/box-drop-large/27000 \
# mittel/2D/table-column/H10 mittel/2D/table-column/H20 mittel/2D/table-column/H40 \
# mittel/3D/table-column/H10 mittel/3D/table-column/H20 mittel/3D/table-column/H40

test_files := $(foreach f, $(tests), $(test_base)/$(f)/$(test_file))
test_xmlfiles := $(foreach f, $(tests), $(test_base)/$(f)/$(test_stem)$(xml_suffix))
test_dirs := $(foreach f, $(tests), $(test_base)/$(f))

sph_octave = $(sph_root)/octave
mfiles     = $(sph_octave)/*.m
mfiles     = 

export SPH_PARAM_B=10000
export SPH_PARAM_BRES=1.5
export SPH_ROOT=$(sph_root)

all: $(test_files)

$(test_base)/2D/pattern-rect-cart/$(test_file): $(scripts)/2D/pattern-rect-cart.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.06 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/pattern-rect-hex/$(test_file): $(scripts)/2D/pattern-rect-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.06 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/pattern-rect-hex-odd/$(test_file): $(scripts)/2D/pattern-rect-hex-odd.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.06 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back


$(test_base)/2D/pattern-test/klein/$(test_file): $(scripts)/2D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.03 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/pattern-test/mittel/$(test_file): $(scripts)/2D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/pattern-test/gross/$(test_file): $(scripts)/2D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/pattern-test/klein/$(test_file): $(scripts)/3D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.03 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/pattern-test/mittel/$(test_file): $(scripts)/3D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/pattern-test/gross/$(test_file): $(scripts)/3D/pattern-test.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

#
# falling-film, monaghan ellipsis(2D), cigar and lentil (3D)
#
$(test_base)/2D/falling-film/klein/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film/mittel/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film/gross/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-hex/klein/$(test_file): $(scripts)/2D/falling-film-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-hex/mittel/$(test_file): $(scripts)/2D/falling-film-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-hex/gross/$(test_file): $(scripts)/2D/falling-film-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-vert/klein/$(test_file): $(scripts)/2D/falling-film-vert.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-vert/mittel/$(test_file): $(scripts)/2D/falling-film-vert.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-vert/gross/$(test_file): $(scripts)/2D/falling-film-vert.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

# falling-film-triang
$(test_base)/2D/falling-film-triang/klein/$(test_file): $(scripts)/2D/falling-film-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-triang/mittel/$(test_file): $(scripts)/2D/falling-film-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/2D/falling-film-triang/gross/$(test_file): $(scripts)/2D/falling-film-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back


$(test_base)/2D/falling-film-2/klein/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/falling-film-2/mittel/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/falling-film-2/gross/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/falling-film-frei/klein/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< none

$(test_base)/2D/falling-film-frei/mittel/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< none

$(test_base)/2D/falling-film-frei/gross/$(test_file): $(scripts)/2D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< none

# falling film 3D
$(test_base)/3D/falling-film/klein/$(test_file): $(scripts)/3D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0125 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/falling-film/mittel/$(test_file): $(scripts)/3D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0075 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/falling-film/gross/$(test_file): $(scripts)/3D/falling-film.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0025 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

# treppe3
$(test_base)/2D/treppe3/klein/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 3 3 0.75

$(test_base)/2D/treppe3/mittel/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 3 3 0.75

$(test_base)/2D/treppe3/gross/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 3 3 0.75

# treppe5
$(test_base)/2D/treppe5/klein/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 5 3 0.75

$(test_base)/2D/treppe5/mittel/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 5 3 0.75

$(test_base)/2D/treppe5/gross/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 5 3 0.75


# treppe11
$(test_base)/2D/treppe11/klein/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 11 6 1.5

$(test_base)/2D/treppe11/mittel/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 11 6 1.5

$(test_base)/2D/treppe11/gross/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 11 6 1.5


# treppe12
$(test_base)/2D/treppe12/klein/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 12 6 1.5

$(test_base)/2D/treppe12/mittel/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 12 6 1.5

$(test_base)/2D/treppe12/gross/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 12 6 1.5


# treppe21
$(test_base)/2D/treppe21/klein/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 21 10 2.5

$(test_base)/2D/treppe21/mittel/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 21 10 2.5

$(test_base)/2D/treppe21/gross/$(test_file): $(scripts)/2D/treppe.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< 21 10 2.5


# treppe3, 3D
$(test_base)/3D/treppe3/klein/$(test_file): $(scripts)/3D/treppe3.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/treppe3/mittel/$(test_file): $(scripts)/3D/treppe3.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/treppe3/gross/$(test_file): $(scripts)/3D/treppe3.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

# treppe5, 3D
$(test_base)/3D/treppe5/klein/$(test_file): $(scripts)/3D/treppe5.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/treppe5/mittel/$(test_file): $(scripts)/3D/treppe5.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back

$(test_base)/3D/treppe5/gross/$(test_file): $(scripts)/3D/treppe5.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< back



$(test_base)/2D/poiseuille-flow/klein/$(test_file): $(scripts)/2D/poiseuille-flow.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/poiseuille-flow/mittel/$(test_file): $(scripts)/2D/poiseuille-flow.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/poiseuille-flow/gross/$(test_file): $(scripts)/2D/poiseuille-flow.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.001 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/poiseuille-flow-hex/klein/$(test_file): $(scripts)/2D/poiseuille-flow-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/poiseuille-flow-hex/mittel/$(test_file): $(scripts)/2D/poiseuille-flow-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both

$(test_base)/2D/poiseuille-flow-hex/gross/$(test_file): $(scripts)/2D/poiseuille-flow-hex.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.001 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< both


$(test_base)/2D/line/mittel/V001/$(test_file): $(scripts)/2D/line.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=1 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/line/mittel/V010/$(test_file): $(scripts)/2D/line.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=10 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/line/mittel/V100/$(test_file): $(scripts)/2D/line.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=100 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/bar/klein/V001/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=1 SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/klein/V010/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=10 SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/klein/V100/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=100 SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/mittel/V001/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=1 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/mittel/V010/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=10 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/mittel/V100/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=100 SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/gross/V001/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=1 SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/gross/V010/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=10 SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/bar/gross/V100/$(test_file): $(scripts)/2D/bar.sh $(mfiles)
	mkdir -p $(@D)
	fluidVel=100 SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/monaghan-ellipsis/klein/$(test_file): $(scripts)/2D/flying-circle.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.03 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/monaghan-ellipsis/mittel/$(test_file): $(scripts)/2D/flying-circle.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/monaghan-ellipsis/gross/$(test_file): $(scripts)/2D/flying-circle.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.005 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-cigar/klein/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.05 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-cigar/mittel/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.025 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-cigar/gross/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-lentil/klein/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.05 SPH_PARAM_INIT_A0_V="100 -100 100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-lentil/mittel/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.025 SPH_PARAM_INIT_A0_V="100 -100 100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/monaghan-lentil/gross/$(test_file): $(scripts)/3D/flying-ball.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_INIT_A0_V="100 -100 100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


#
# table-shot 2D and 3D
#
$(test_base)/2D/table-shot/klein/V001/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) fluidVel="1" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/klein/V010/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) fluidVel="10" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/klein/V100/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) fluidVel="100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/mittel/V001/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) fluidVel="1" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/mittel/V010/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) fluidVel="10" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/mittel/V100/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) fluidVel="100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/gross/V001/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) fluidVel="1" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/gross/V010/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) fluidVel="10" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-shot/gross/V100/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) fluidVel="100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/3D/table-shot/klein/V001/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 fluidVel="1" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-shot/klein/V010/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 fluidVel="10" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-shot/klein/V100/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 fluidVel="100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-shot/mittel/V001/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 fluidVel="1" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-shot/mittel/V010/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 fluidVel="10" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-shot/mittel/V100/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 fluidVel="100" SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<



$(test_base)/2D/table-column/klein/H10/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.10 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-column/klein/H20/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.20 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-column/klein/H40/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.40 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-column/mittel/H10/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.10 SPH_PARAM_DX=0.002 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-column/mittel/H20/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.20 SPH_PARAM_DX=0.002 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/table-column/mittel/H40/$(test_file): $(scripts)/2D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.40 SPH_PARAM_DX=0.002 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/3D/table-column/klein/H10/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.10 SPH_PARAM_DX=0.05 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-column/klein/H20/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.20 SPH_PARAM_DX=0.05 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-column/klein/H40/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.40 SPH_PARAM_DX=0.05 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-column/mittel/H10/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.10 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-column/mittel/H20/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.20 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-column/mittel/H40/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	CHEIGHT=0.40 SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/box-drop/klein/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-drop/mittel/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-drop/gross/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/box/klein/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< no

$(test_base)/2D/box/mittel/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< no

$(test_base)/2D/box/gross/$(test_file): $(scripts)/2D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $< no


$(test_base)/2D/utube/klein/$(test_file): $(scripts)/2D/utube.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/utube/mittel/$(test_file): $(scripts)/2D/utube.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/utube/gross/$(test_file): $(scripts)/2D/utube.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/artesian-fountain/klein/$(test_file): $(scripts)/2D/artesian-fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/artesian-fountain/mittel/$(test_file): $(scripts)/2D/artesian-fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/artesian-fountain/gross/$(test_file): $(scripts)/2D/artesian-fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/3D/box-drop/klein/$(test_file): $(scripts)/3D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/box-drop/mittel/$(test_file): $(scripts)/3D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.02 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/box-drop/gross/$(test_file): $(scripts)/3D/box-drop.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/fountain/klein/$(test_file): $(scripts)/2D/fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/fountain/mittel/$(test_file): $(scripts)/2D/fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/fountain/gross/$(test_file): $(scripts)/2D/fountain.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/fountain/klein/$(test_file): $(scripts)/3D/fountain2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.03 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/fountain/mittel/$(test_file): $(scripts)/3D/fountain2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.015 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/fountain/gross/$(test_file): $(scripts)/3D/fountain2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0075 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/dambreak/klein/$(test_file): $(scripts)/2D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak/mittel/$(test_file): $(scripts)/2D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak/gross/$(test_file): $(scripts)/2D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak/riesig/$(test_file): $(scripts)/2D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxRiesig2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/3D/dambreak/klein/216/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.081 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/klein/1000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.056 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/klein/2744/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/mittel/10000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.027 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/mittel/50000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.016 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/gross/200000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/gross/1000000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.006 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak/gross/2000000/$(test_file): $(scripts)/3D/dambreak.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0047 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/2D/dambreak-obstacle/klein/$(test_file): $(scripts)/2D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.08 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak-obstacle/mittel/$(test_file): $(scripts)/2D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak-obstacle/gross/$(test_file): $(scripts)/2D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.02 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/dambreak-obstacle/riesig/$(test_file): $(scripts)/2D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak-obstacle/klein/$(test_file): $(scripts)/3D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.08 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak-obstacle/mittel/$(test_file): $(scripts)/3D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak-obstacle/gross/$(test_file): $(scripts)/3D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.02 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/dambreak-obstacle/riesig/$(test_file): $(scripts)/3D/dambreak-obstacle-triang.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/klein/3D/table/27/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0334 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/klein/3D/table/125/$(test_file): $(scripts)/3D/szene-tafel2.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.019 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/3D/table/10000/$(test_file): $(scripts)/3D/szene-tafel1.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/3D/table/50000/$(test_file): $(scripts)/3D/szene-tafel1.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.005 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<



$(test_base)/3D/table-skewed/klein/27/$(test_file): $(scripts)/3D/szene-tafel-schraeg.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.003 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-skewed/klein/125/$(test_file): $(scripts)/3D/szene-tafel-schraeg.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0022 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-skewed/mittel/10000/$(test_file): $(scripts)/3D/szene-tafel-schraeg.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0005 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/3D/table-skewed/mittel/50000/$(test_file): $(scripts)/3D/szene-tafel-schraeg.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0005 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<



$(test_base)/2D/box-filled/winzig/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.015 wallOffset=0 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled/klein/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.008 wallOffset=0 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled/mittel/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0038 wallOffset=0 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled/gross/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0022 wallOffset=0 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled-offset/klein/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.008 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled-offset/mittel/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0038 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/2D/box-filled-offset/gross/$(test_file): $(scripts)/2D/box-filled.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.0022 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/mittel/2D/box-filled-large/121/$(test_file): $(scripts)/2D/box-filled-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.08 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-filled-large/625/$(test_file): $(scripts)/2D/box-filled-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-filled-large/2000/$(test_file): $(scripts)/2D/box-filled-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.02 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-filled-large/10000/$(test_file): $(scripts)/2D/box-filled-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-drop-large/121/$(test_file): $(scripts)/2D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.08 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-drop-large/625/$(test_file): $(scripts)/2D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.04 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-drop-large/2000/$(test_file): $(scripts)/2D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.02 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/2D/box-drop-large/10000/$(test_file): $(scripts)/2D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.01 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/mittel/3D/box-drop-large/261/$(test_file): $(scripts)/3D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.15 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/3D/box-drop-large/1728/$(test_file): $(scripts)/3D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.08 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/3D/box-drop-large/8000/$(test_file): $(scripts)/3D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.048 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/mittel/3D/box-drop-large/27000/$(test_file): $(scripts)/3D/box-drop-large.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=0.032 SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


#
$(test_base)/1D/burgers/klein/$(test_file): $(scripts)/1D/burgers.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<


$(test_base)/1D/stempel/klein/$(test_file): $(scripts)/1D/stempel.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxKlein2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/1D/stempel/mittel/$(test_file): $(scripts)/1D/stempel.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxMittel2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

$(test_base)/1D/stempel/gross/$(test_file): $(scripts)/1D/stempel.sh $(mfiles)
	mkdir -p $(@D)
	SPH_PARAM_DX=$(dxGross2) SPH_PARAM_SCENE_OUTFILE=$(basename $@) $<

PREFIX ?= $(HOME)/local
INSTDIR = $(PREFIX)/sph-data
INSTALL ?= install
CP ?= cp
MKDIR_P ?= mkdir -p

install:
	$(MKDIR_P) $(INSTDIR)/test-cases2
	true $(foreach f, $(test_files), && ( \
		echo "$(f) -> $(INSTDIR)/$(dir $(f))"; \
		echo "$(f).env -> $(INSTDIR)/$(dir $(f))"; \
		echo "$(basename $(f)).sphscene.xml -> $(INSTDIR)/$(dir $(f))"; \
		$(MKDIR_P) $(INSTDIR)/$(dir $(f)); \
		$(CP) $(f) $(f).env $(basename $(f)).sphscene.xml $(INSTDIR)/$(dir $(f)); ))

clean:
	$(RM)  $(wildcard core* *.vtk *.env boundary*.txt *.sphscene.xml)

clean-tests:
	$(RM) -r $(test_dirs)

check: $(test_xmlfiles)
	@for i in $(test_xmlfiles); do \
		echo $$i; \
		xsltproc -o /dev/null ../xsl/scene-check.xsl $$i; \
	done

info: $(test_xmlfiles)
	@for i in $(test_xmlfiles); do \
		echo $$i; \
		xsltproc ../xsl/scene-info.xsl $$i; \
	done

BUNZIP2 ?= bunzip2
TAR ?= tar
unpack:
	$(BUNZIP2) -cd test-cases.tar.bz2 | $(TAR) -xf -
	$(BUNZIP2) -cd test-cases-contrib.tar.bz2 | $(TAR) -xf -
