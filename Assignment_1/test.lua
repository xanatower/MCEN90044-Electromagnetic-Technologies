showconsole()
mydir="./"
open(mydir .. "whole.fem")
mi_saveas(mydir .. "temp1.fem")
mi_seteditmode("group")
for n=0,15 do
    mi_analyze()
    mi_loadsolution()
    mo_groupselectblock(1)
    fz=mo_blockintegral(19)
    print(fz)
    if (n<15) then
        mi_selectgroup(3)
        mi_movetranslate(0,0.1)
        a
    end
end
mo_close()
mi_close()


print(clock());
if clock_update == false then
    clock_now =  clock();
    clock_diff = clock_now - clock_before;
    