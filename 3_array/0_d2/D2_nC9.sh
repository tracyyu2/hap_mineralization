20250403
/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-10 sym=D2_Z.sym -out:suffix _D2neg10 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all


/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/fiber_design/Dihedral/00_dihedral_rosetta/input/d5_D9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-3.0 sym=D3_Z.sym angle=0 nres=795 -out:suffix _trans-3.0_rot0 -parser:protocol dihedral_symmetrize_and_translate_and_spin.xml -mute all


/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-1 sym=D2_Z.sym -out:suffix _D2neg1 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-30 sym=D2_Z.sym -out:suffix _D2neg30 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-35 sym=D2_Z.sym -out:suffix _D2neg35 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-40 sym=D2_Z.sym -out:suffix _D2neg40 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/d3_nC9_chainA_Ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-45 sym=D2_Z.sym -out:suffix _D2neg45 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

20250404

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/nC9_monomer_d2_aligned.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-5 sym=D2_Z.sym -out:suffix _D2neg5 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all

/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease -database /software/rosetta/latest/database @rosetta_scripts.flag -s /home/tracyy2/array_design/00_d2_rosetta_symmetrize/input/nC9_monomer_d2_aligned_ala.pdb -scorefile score.sc  -parser:script_vars translate_X=0 translate_Y=0 translate_Z=1 distance=-6.5 sym=D2_Z.sym -out:suffix _D2neg1 -parser:protocol dihedral_symmetrize_and_translate.xml -mute all


