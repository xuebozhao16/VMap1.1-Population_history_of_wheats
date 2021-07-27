import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

mu = 1.25e-8
gen = 1
#dir_ = "MSMC2_OUTPUT"
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/wild_emmer.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu


#plt.figure(figsize=(8, 10))
#plt.subplot(211)
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='red', label='Yoruba')
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='French')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
#plt.subplot(212)
relativeCCR=2.0 * msmc_out['lambda'] / (msmc_out['lambda'] + msmc_out['lambda'])
plt.semilogx(t_years,relativeCCR, drawstyle='steps')
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()
#plt.savefig("MSMC_plot.pdf")

#####现在做的是四倍体整体的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Free threshing')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Wild emmer')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Dom emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Tetraploid_1.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Tetraploid_2.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_11 + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Tetraploid_3.pdf')
plt.show()

#####现在做的是六倍体AB整体的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine2_AB.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceAB')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EU landraceAB')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EA landraceAB')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Hexaploid_1.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Hexaploid_2.pdf')
plt.show()

###现在想把四倍体和六倍体合起来
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Free threshing')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Wild emmer')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Dom emmer')
msmc_out2=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine2_AB.txt", sep='\t', header=0)
t_years2=gen * ((msmc_out2.left_time_boundary + msmc_out2.right_time_boundary)/2) / mu
plt.semilogx(t_years2, (1/msmc_out2.lambda_00)/(2*mu), drawstyle='steps',color='green', label='WA landraceAB')
plt.semilogx(t_years2, (1/msmc_out2.lambda_01)/(2*mu), drawstyle='steps',color='pink', label='EU landraceAB')
plt.semilogx(t_years2, (1/msmc_out2.lambda_11)/(2*mu), drawstyle='steps',color='yellow', label='EA landraceAB')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()



#####现在做的是D的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine2_D2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EU landraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
#relativeCCR00=2.0 * msmc_out.lambda_00 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
#relativeCCR11=2.0 * msmc_out.lambda_11 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
#plt.semilogx(t_years,relativeCCR00, color='red',drawstyle='steps')
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
#plt.semilogx(t_years,relativeCCR11, color='Orange',drawstyle='steps')
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/diploid_3.pdf')
plt.show()

#####现在做的是D的情况，加上strangulta
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/combine2_D.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EU landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Strangulata')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Hexaploid_3.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Hexaploid_4.pdf')
plt.show()

def getCCRintersect(df, val):
    xVec = gen * ((df.left_time_boundary + df.right_time_boundary)/2) / mu
    yVec = 2.0 * df.lambda_01 / (df.lambda_00 + df.lambda_11)
    i = 0
    while yVec[i] < val:
        i += 1
    assert i > 0 and i <= len(yVec), "CCR intersection index out of bounds: {}".format(i)
    assert yVec[i - 1] < val and yVec[i] >= val, "this should never happen"
    intersectDistance = (val - yVec[i - 1]) / (yVec[i] - yVec[i - 1])
    return xVec[i - 1] + intersectDistance * (xVec[i] - xVec[i - 1])


print(getCCRintersect(msmc_out, 0.5)) #Print out the time when relativeCCR=0.5


mu = 6.5e-9
gen = 1

##Urartu
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8/wild_emmer.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Urartu')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
##wild emmer
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/wild_emmer.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='green', label='Wild emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
##dom emmer
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8/dom_emmer.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Wild emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
##dom emmer
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8/Rivet_wheat.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Wild emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
####现在是8个haplotype的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8/MSMC_AB_1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EU landraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
#relativeCCR00=2.0 * msmc_out.lambda_00 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
#relativeCCR11=2.0 * msmc_out.lambda_11 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
#plt.semilogx(t_years,relativeCCR00, color='red',drawstyle='steps')
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
#plt.semilogx(t_years,relativeCCR11, color='Orange',drawstyle='steps')
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()


###继续画图,现在做的的haplotype8_plot2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

mu = 6.5e-9
gen = 1
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/wild_emmer.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu

##二倍体
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_3.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Urartu')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Speltoids')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Wild emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/diploid_1.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()
#####现在做的是四倍体整体的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_AB_1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Rivet wheat')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Wild emmer')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Dom emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()

#####现在做的是四倍体整体的情况2
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_AB_2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Persion wheat')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Rivet wheat')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Durum')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()

#####现在做的是六倍体AB整体的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_AB_1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='LandraceAB')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Indian_DwarfAB')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='YunanAB')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/Tetraploid_2.pdf')
plt.show()

#####现在做的是六倍体AB整体的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_AB_2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceAB')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EU landraceAB')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EA landraceAB')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V1/diploid_2.pdf')
plt.show()

#####现在做的是D的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_ABD_D2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='YunanD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Indian_DwarfABD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='LandraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()

#####现在做的是D的情况
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/MSMC_ABD_D1.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EU landraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11  + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.show()

##Urartu
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/Urartu_2.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Urartu')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()

##Speltoides
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/speltoides_3.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Speltoides')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()

##Strangulata
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/hap8_V2/Strangulata_3.final.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out['lambda'])/(2*mu), drawstyle='steps',color='blue', label='Strangulata')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()


####现在做没有Introgression的情况 AB
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/MSMC_AB3.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EU landraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_01 + msmc_out.lambda_11)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/WE_North_DE.pdf')
plt.show()




####现在做没有Introgression的情况 ABD
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/MSMC_AB_ABD2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='WA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='EA landraceD')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='EU landraceD')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_01 + msmc_out.lambda_11)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/DE_Free.pdf')
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_11 / (msmc_out.lambda_00 + msmc_out.lambda_01 + msmc_out.lambda_11)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/Free_breadwheat.pdf')
plt.show()

##二倍体
msmc_out=pd.read_csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/MSMC_A2.txt", sep='\t', header=0)
t_years=gen * ((msmc_out.left_time_boundary + msmc_out.right_time_boundary)/2) / mu
plt.semilogx(t_years, (1/msmc_out.lambda_00)/(2*mu), drawstyle='steps',color='red', label='Urartu')
plt.semilogx(t_years, (1/msmc_out.lambda_01)/(2*mu), drawstyle='steps',color='blue', label='Speltoids')
plt.semilogx(t_years, (1/msmc_out.lambda_11)/(2*mu), drawstyle='steps',color='Orange', label='Wild emmer')
plt.xlabel("years ago")
plt.ylabel("population Sizes")
plt.legend()
plt.show()
relativeCCR01=2.0 * msmc_out.lambda_01 / (msmc_out.lambda_00 + msmc_out.lambda_11 + msmc_out.lambda_01)
plt.semilogx(t_years,relativeCCR01, color='blue',drawstyle='steps')
plt.hlines(0.5, 0, 1000000, colors = "c", linestyles = "dashed")
plt.xlabel("years ago")
plt.ylabel("Relative CCR")
plt.savefig('/Users/xuebozhao/Documents/LuLab/wheatSpeciation/MSMC/V2_no_introgression/WEI_DEI.pdf')
plt.show()





