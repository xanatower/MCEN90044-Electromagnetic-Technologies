-- Example 7.1 from Bauer
-- gregheins@ieee.org

-- Note: Requires FEMM 4.2

-- set up model
steps = 8;
stepIncrement = 0.5;
initialLocation = 0.0;

-- load geometry
mydir="./"
open(mydir .. "Afternoon_exercise.FEM")

-- and save as a temporary file name so we don't 
-- ovewrite the original geometry
mi_saveas(mydir .. "Afternoon_temp.fem")

-- show the Lua console window so that we can see the 
-- program report its progress
showconsole()

-- move the motion region
-- save parameters of interest at each point

d = {};

for k = 1,(steps+1) do
	
  	print((k-1) .. "/" .. steps);
  	d[k]={};
  	d[k][1] = initialLocation+(k-1)*stepIncrement;

  	-- compute the flux linkage at three different currents
  	
 	mi_analyze(1);
  	mi_loadsolution();

 	--d[k][2],v2,d[k][3] = mo_getcircuitproperties("Coil"); -- current, voltage, flux linkage;
 	A, d[k][2], d[k][3], Sig, E, H1, H2, Je, Js, Mu1, Mu2, Pe, Ph =  mo_getpointvalues(2, -3); -- current, voltage, flux linkage
  	-- increment the motion
 	mi_selectgroup(1);
	mi_movetranslate(0, stepIncrement, 4); -- 4 = group (p 87 of FEMM manual)
  	mi_clearselected();
end


-- write results to disk for further analysis.
fp=openfile(mydir .. "Example 11-3.csv","w")
for k = 1,(steps+1) do
	--d[k][1],",",d[k][2],",",d[k][3]-> displacement, Bx, By
  write(fp,d[k][1],",",d[k][2],",",d[k][3],",","\n")
end
closefile(fp);