.PHONY: run build show

build:
	scons

run:
	cd output; ../build/ADV1D

show:
	paraview --state=plot.pvsm
	# paraview --data=output/adv1d_..vtr
