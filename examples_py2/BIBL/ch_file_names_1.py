import os
os.chdir("c:/bibl/pdf/cur_op_struct_biol")
fnames = os.listdir(".")
for fn in (fnames):
  fn2 = fn
  fn2 = fn2.replace("cur_op","curr_opin")
  print(fn2)
  os.rename(fn,fn2)


