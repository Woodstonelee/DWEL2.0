function dt2nb, data_type
compile_opt idl2
;
;
numbyt=-1
case data_type of
0:numbyt=0
1:numbyt=1
2:numbyt=2
3:numbyt=4
4:numbyt=4
5:numbyt=8
6:numbyt=8
9:numbyt=16
12:numbyt=2
13:numbyt=4
14:numbyt=8
15:numbyt=8
else:numbyt=-1
endcase
;
return,numbyt
end
