--------------------------------------------
-- outrunner motor example
-- greg@heins.com.au, 20187
-- based on script from edie_currants@yahoo.com  
-- Note: Requires FEMM 4.2
--------------------------------------------

-- set up model
steps_per_electrical_cycle = 18; --should be multiple of 3
number_of_electrical_cycles = 1;

number_of_poles = 14; 
phase_current_rms = 0; 

steps = steps_per_electrical_cycle*number_of_electrical_cycles;
current_mag = sqrt(2)*phase_current_rms;
current_offset_angle = -120; --aline current with back emf

-- load geometry
mydir="./"
open(mydir .. "untitled_part2_90.FEM")

-- and save as a temporary file name so we don't 
-- ovewrite the original geometry
mi_saveas(mydir .. "temp.fem")

-- show the Lua console window so that we can see the 
-- program report its progress
showconsole()

-- move the rotor through one pole pitch in small
-- increments, recording the flux linkage of each
-- phase at each rotor position
d = {};-- to store things

for k = 1,(steps+1) do
	angle_elec_degrees_increment = 360/steps_per_electrical_cycle;
	angle_mech_degrees_increment = angle_elec_degrees_increment/(number_of_poles/2);
	angle_elec_degrees = (k-1)*angle_elec_degrees_increment;
	angle_mech_degrees = angle_elec_degrees/(number_of_poles/2);

	
  	mi_modifycircprop("A",1,current_mag*cos(rad(angle_elec_degrees+current_offset_angle-120)));
  	mi_modifycircprop("B",1,current_mag*cos(rad(angle_elec_degrees+current_offset_angle)));
  	mi_modifycircprop("C",1,current_mag*cos(rad(angle_elec_degrees+current_offset_angle+120)));
	
  	print((k-1) .. "/" .. steps);
  	mi_analyze(1);
  	mi_loadsolution();
  	d[k]={};
  	d[k][1] = angle_elec_degrees;

  	-- collect the flux linkage and current for each phase
   	d[k][6],v2,d[k][3] = mo_getcircuitproperties("A"); -- current, voltage, flux linkage
    d[k][7],v2,d[k][4] = mo_getcircuitproperties("B"); 
    d[k][8],v2,d[k][5] = mo_getcircuitproperties("C");

  	-- compute the torque
 	mo_groupselectblock(2);
 	d[k][2] = mo_blockintegral(22);

  	-- increment the rotor's position
  	mi_selectgroup(2);
  	mi_moverotate(0, 0, angle_mech_degrees_increment , 4);
  	mi_clearselected();
end

-- write results to disk for further analysis.
fp=openfile(mydir .. "Cogging_Torque.csv_3","w")
for k = 1,(steps+1) do
   write(fp,d[k][1],",",d[k][2],",",d[k][3],",",d[k][4],",",d[k][5],",",d[k][6],",",d[k][7],",",d[k][8],"\n")
end
closefile(fp);
