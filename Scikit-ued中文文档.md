This document is translated from

> http://scikit-ued.readthedocs.io/en/master/index.html

THANK YOU! 



<!-- TOC -->

- [Tutorials](#tutorials)
    - [1 建立原子模型](#1-建立原子模型)
        - [1.1 The Crystal Class](#11-the-crystal-class)
            - [1.1.1 构造Crystal](#111-构造crystal)
            - [1.1.2 手动构造Crystal](#112-手动构造crystal)
            - [1.1.3 Crystal attributes](#113-crystal-attributes)
            - [1.1.4 Lattice vectors and reciprocal space](#114-lattice-vectors-and-reciprocal-space)
            - [1.1.5 Space-group Information](#115-space-group-information)
            - [1.1.6 和ASE的兼容](#116-和ase的兼容)
        - [1.2 The Atom Class](#12-the-atom-class)
    - [2 Simulation Tutorial](#2-simulation-tutorial)
        - [2.1 Polycrystalline 多晶衍射模拟](#21-polycrystalline-多晶衍射模拟)
        - [2.2 Electrostatic 静电电位模拟](#22-electrostatic-静电电位模拟)
    - [3 Baseline Tutorials](#3-baseline-tutorials)
        - [3.1 使用离散小波变换的迭代基线确定](#31-使用离散小波变换的迭代基线确定)
        - [3.2 使用双树复小波变换的迭代基线确定](#32-使用双树复小波变换的迭代基线确定)
    - [4 图像分析与处理](#4-图像分析与处理)
        - [4.1 衍射图像对齐](#41-衍射图像对齐)
        - [4.2 涉及对称性的图像处理](#42-涉及对称性的图像处理)
        - [4.3 Pixel Masks](#43-pixel-masks)
            - [4.3.1 创建一个pixel mask](#431-创建一个pixel-mask)
        - [4.4 多晶衍射图的图像分析](#44-多晶衍射图的图像分析)
            - [4.4.1 确定中心](#441-确定中心)
            - [4.4.2 角度分布](#442-角度分布)
        - [4.5 Bonus:去除热点](#45-bonus去除热点)
    - [5 Plotting Tutorial](#5-plotting-tutorial)
        - [color](#color)
- [Reference/API](#referenceapi)
    - [1 Baseline-determination](#1-baseline-determination)
        - [1.1 skued.baseline_dt](#11-skuedbaseline_dt)

<!-- /TOC -->



#  Tutorials

Tutorials用于基本操作，如果需要更具体的操作，请查看 Reference/API

## 1 建立原子模型

### 1.1 The Crystal Class

依赖于 `skued.structure`包

#### 1.1.1 构造Crystal

```py
from skued import Crystal
TiSe2 = Crystal.from_cif('tise2.cif')
# Crystal Information File (CIF, .cif)
```

Scikit-ued同样拥有外部数据库的CIF文件，有效名在`Crystal.builtins`中

```py
assert 'Au' in Crystal.builtins
gold = Crystal.from_database('Au')
```

对于Protein DataBank文件，只需要提供4个字母的身份code：

```py
hemoglobin = Crystal.from_pdb('1gzx')
```

也可通过`Crystallography Open Database`构造`Crystal`

```py
# Default is the latest revision
vo2 = Crystal.from_cod(1521124)

# Revisions are accessible as well
old_vo2 = Crystal.from_cod(1521124, revision = 140771)
```

#### 1.1.2 手动构造Crystal

如果没有文件，或者想要构造一个理想化晶体，可以手动构造。

构造需要
1. iterable of Atom objects--必须是满单元结构
2. three lattice vectors

创建一个最简单的晶体结构： `alpha-Polonium (simple cubic)`：
```py
from skued import Crystal, Atom
import numpy as np

lattice_vectors = 3.35 * np.eye(3)
unitcell = [Atom('Po', coords = [0,0,0])]

polonium = Crystal(unitcell, lattice_vectors)
```

#### 1.1.3 Crystal attributes

首先，`Crystal`是个iterable

```py
from skued import Crystal
graphite = Crystal.from_database('C')

for atm in graphite:    #Loops over atoms in the unit cell
    print(atm)
```

The len() of a Crystal is the unit cell size (in number of atoms)：
```py
c = Crystal.from_pdb('1gzx') # hemoglobin
len(c)                          # 17536
```

`Crystal`是set-like容器，检查容器关系(使用内置in申明)：
```py
graphite = Crystal.from_database('C')
carbon = next(iter(graphite))

assert carbon in graphite
```

`Crystal` instances can be equated to each other:
```py
gold = Crystal.from_database('Au')
silver = Crystal.from_database('Ag')

assert gold == silver # false
```

如果`Crystal` 是从文件中产生的，文件路径可以从`source`中获取：
```py
c = Crystal.from_pdb('1gzx')
print(c.source)
```

`Crystal` 实例有一个很好的字符串表达式：
```py
lsmo = Crystal.from_database(‘LSMO’) print(lsmo)
```

`Crystal` 也可被转化为`NumPy`矩阵
```py
import numpy
arr = numpy.array(Crystal.from_database('Si'))
    print(arr)
```
#### 1.1.4 Lattice vectors and reciprocal space

[What is the difference between lattice vectors and basis vectors (in crystallography)?](https://www.quora.com/What-is-the-difference-between-lattice-vectors-and-basis-vectors-in-crystallography)


当创建完`Crystal`后，可以使用隐含的Lattice superclass来控制lattice参数


```py
# 内置的石墨(graphite)示例
from skued import Crystal
graphite = Crystal.from_database('C')

a1, a2, a3 = graphite.lattice_vectors
b1, b2, b3 = graphite.reciprocal_vectors

# 标准的三长度和三角度用于描述lattice
a, b, c, alpha, beta, gamma = graphite.lattice_parameters

# unit cell volume
vol = graphite.volume   # Angstroms cubed
density = vol/len(graphite)
```

#### 1.1.5 Space-group Information

```py
# 可以通过lattice_system属性来获取the lattice system of a Lattice or Crystal instance

vo2 = Crystal.from_database('vo2-m1') # Monoclinic M1 VO2  单斜(晶系)的
print(vo2.lattice_system)             # = 'monoclinic'

# 可以通过lattice_system()来对length tolerances（这是什么？？）进行更好的控制
# 可以从Crystal实例中获取 space-group information

from skued import Crystal

gold = Crystal.from_database('Au')
spg_info = gold.spacegroup_info()

# 在上述的例子中，spg_info是有下述key的目录

# 'international_symbol': International Tables of Crystallography space-group symbol (short);
# 'international_full': International Tables of Crystallography space-group full symbol;
# 'hall_symbol' : Hall symbol;
# 'pointgroup' : International Tables of Crystallography point-group;
# 'international_number' : International Tables of Crystallography space-group number (between 1 and 230);
# 'hall_number' : Hall number (between 1 and 531).

# Scattering utilities散射实用程序
# 当处理散射数据和建模时，lattice模型有一些简化方法

# Miller indices和 scattering vectors 的转化：
from skued import Crystal
graphite = Crystal.from_database('C')

# Behavior inherited from Lattice superclass
G = graphite.scattering_vector(1,0,0)
h, k, l = graphite.miller_indices(G) #1, 0, 0

```

#### 1.1.6 和ASE的兼容


### 1.2 The Atom Class



```py
# To create an atom, simply provide its element and coordinates:

from skued import Atom
copper = Atom(element = 'Cu', coords = [0,0,0])



from skued import Crystal
graphite = Crystal.from_database('C')
carbon = list(graphite)[-1]
fractional = carbon.coords
real = carbon.xyz(lattice = graphite)




copper = Atom('Cu', coords = [0,0,0])
silver = Atom('Ag', coords = [1,0,0])
dist = silver - copper                  # distance in fractional coordinates
```




## 2 Simulation Tutorial

模拟衍射部分。

### 2.1 Polycrystalline 多晶衍射模拟 

准备内容：一个`Crystal`和一个散射长度范围，s as s=sinθ/λ:

```py

import matplotlib.pyplot as plt
import numpy as np
from skued import Crystal
graphite = Crystal.from_database('C')
from skued import powdersim
s = np.linspace(0.1, 0.8, 1000)
diff = powdersim(crystal = graphite, scattering_length = s)
plt.figure()
plt.plot(s, diff/diff.max())
plt.xlim([s.min(), s.max()])
plt.xlabel('$s = \sin{\\theta}/\lambda$')
plt.ylabel('Diffracted intensity (A.u.)')
plt.title('Polycrystalline graphite diffraction')

```

### 2.2 Electrostatic 静电电位模拟


电子的散射电位是晶体的静电势，因此计算这个电势是很有用的。

```py
# electrostatic()用于无限空间
from skued import Crystal
from skued import electrostatic

graphite = Crystal.from_database('C')

extent = np.linspace(-10, 10, 256)
xx, yy, zz = np.meshgrid(extent, extent, extent)
potential = electrostatic(graphite, xx, yy, zz)
```

```py
# pelectrostatic()用于x-y无限，z有限

from skued import pelectrostatic

extent = np.linspace(-5, 5, 256)
xx, yy= np.meshgrid(extent, extent)
potential = pelectrostatic(graphite, xx, yy)

```



## 3 Baseline Tutorials

电子衍射的背景信号比X射线衍射要多很多。

在多晶样品（即1D衍射信号）的情况下，去除背景信号的确定方式是使用基于双树复数小波变换的迭代方法。

```py
import matplotlib.pyplot as plt
import numpy as np

s, intensity = np.load('powder.npy')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(s, intensity, 'k')
ax.set_xlabel('s ($4 \pi / \AA$)')
ax.set_ylabel('Diffracted intensity (counts)')
ax.set_title('Background-subtracted diffraction pattern of rutile VO$_2$')
plt.show()
```


从衍射图案的无弹性散射区域内插背景将非常困难，而且这是适度的我们可以添加典型的氮化硅衬底背景以及非弹性散射效应：

```py
from skued import gaussian

background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

```
### 3.1 使用离散小波变换的迭代基线确定

> Galloway et al. ‘An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet Transforms’, Applied Spectroscopy pp. 1370 - 1376, September 2009.




```py
import numpy as np
from skued import gaussian
from skued import baseline_dwt

s, intensity = np.load('powder.npy')

# Double exponential inelastic background and substrate effects
diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

baseline = baseline_dwt(signal, level = 6, max_iter = 150, wavelet = 'sym6')
```


### 3.2 使用双树复小波变换的迭代基线确定

对于 1D数据，双树复小波变换比`baseline_dwt()`性能要高很多。

> L. P. René de Cotret and B. J. Siwick, A general method for baseline-removal in ultrafast electron powder diffraction data using the dual-tree complex wavelet transform, Struct. Dyn. 4 (2017)

```py
import numpy as np
from skued import gaussian
from skued import baseline_dt

s, intensity = np.load('powder.npy')

# Double exponential inelastic background and substrate effects
diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

baseline = baseline_dt(signal, wavelet = 'qhisft3', level = 6, max_iter = 150))
```

`baseline_dt()`比`baseline_dwt()`更精准，但是`baseline_dwt()`可以用于1D和2D的数据。



## 4 图像分析与处理

衍射图分析本质上是专门的图像处理。主要依赖于`skued.image`+`npstreams`

### 4.1 衍射图像对齐

衍射图样可能是几分钟内的，为了可靠的数据合成，基于参考线对齐很重要

通常通过测量图像之间的cross-correlation来检测、登记两个相似图像之间的转换。 当图像非常相似时，只需要看看scikit-image的skimage.feature.register_translation（）。

然而，衍射图样有一个固定特征：beam-block的位置。 因此，在计算cross-correlation时，必须忽略衍射图案中的一些像素。

将“无效像素”设置为0将不起作用。 我们必须通过scikit-ued的mnxc2（）来使用masked normalized cross-correlation.

所有这些都在scikit-ued的diff_register（）函数中处理。 我们来看看一些多晶铬：



```py
from skued import diffread
import matplotlib.pyplot as plt

ref = diffread('Cr_1.tif')
im = diffread('Cr_2.tif')

fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize = (9,3))
ax1.imshow(ref, vmin = 0, vmax = 200)
ax2.imshow(im, vmin = 0, vmax = 200)
ax3.imshow(ref - im, cmap = 'RdBu_r', vmin = -100, vmax = 100)

for ax in (ax1, ax2, ax3):
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

ax1.set_title('Reference')
ax2.set_title('Data')
ax3.set_title('Difference')

plt.tight_layout()
plt.show()
```

为了决定额外的相移，使用mask来使beam-block 和main-beam模糊。


```py
from skued import diff_register, shift_image, diffread
import numpy as np

ref = diffread('Cr_1.tif')
im = diffread('Cr_2.tif')

mask = np.zeros_like(ref, dtype = np.bool)
mask[0:1250, 950:1250] = True

shift = diff_register(im, reference = ref, mask = mask)
im = shift_image(im, shift)
```



### 4.2 涉及对称性的图像处理


基于晶体结构衍射图呈现旋转对称性。 我们可以利用这种对称性在出现伪像或缺陷的情况下校正图像。 一个有用的例程是nfold（），它基于旋转对称对衍射图案的各个部分进行平均。

使用nfold()，我们需要知道衍射图案的中心。

```py
from skued import nfold, diffread

im = diffread('graphite.tif')
av = nfold(im, mod = 6, center = center)    # mask is optional
```


```py
# 完整版代码
import matplotlib.pyplot as plt
from skued import nfold, diffread
import numpy as np

center = (1010, 1111)

mask = np.zeros((2048, 2048), dtype = np.bool)
mask[1100::, 442:480] = True # Artifact line
mask[0:1260, 900:1140] = True # beamblock

image = diffread('graphite.tif')
av = nfold(image, mod = 6, center = center, mask = mask)

fig , (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (9,3))
ax1.imshow(image, vmin = 0, vmax = 150)
ax2.imshow(mask, vmin = 0, vmax = 1)
ax3.imshow(av, vmin = 0, vmax = 150)

for ax in (ax1, ax2, ax3):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

ax1.set_title('Graphite')
ax2.set_title('Mask')
ax3.set_title('Averaged')

plt.tight_layout()
plt.show()
```



### 4.3 Pixel Masks

pixel masks 可以拒绝图像数据的单个像素。这些masks由无效像素的布尔矩阵评估为True来表示。

#### 4.3.1 创建一个pixel mask

使用`dark_run_*.tif.`文件来展示这一系列的图像。通过`mask_from_collection()`来创建一个pixel mask.

```py
from glob import iglob
from skued import mask_from_collection, diffread

dark_runs = map(diffread, iglob('dark_run_*.tif'))    # Can be a huge stack of images
mask = mask_from_collection(dark_runs)
```

在上例中，[0,30000]范围**内?(不确定)**的像素值将被标记为无效（默认行为）。 此外，计算图像集上的每像素标准差; 波动太大的像素也被拒绝。


### 4.4 多晶衍射图的图像分析

#### 4.4.1 确定中心

多晶衍射图显示同心环，找到那些同心环的中心是重要的。

```py
from skued import diffread
import matplotlib.pyplot as plt
path = 'Cr_1.tif'

im = diffread(path)
mask = np.zeros_like(im, dtype = np.bool)
mask[0:1250, 950:1250] = True

im[mask] = 0
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(im, vmin = 0, vmax = 200)

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()
```
这个多晶二氧化钒的衍射图案非常嘈杂。 使用`powder_center()`可以找到这种对称模式的中心

```py
from skued import powder_center, diffread
import numpy as np
import matplotlib.pyplot as plt
path = 'Cr_1.tif'
im = diffread(path)
mask = np.zeros_like(im, dtype = np.bool)
mask[0:1250, 950:1250] = True
ic, jc = powder_center(im, mask = mask)
ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),indexing = 'ij')
rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
im[rr < 100] = 1e6

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(im, vmin = 0, vmax = 200)

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

plt.show()
```

#### 4.4.2 角度分布

首先，创建一个测试图像：


```py
import numpy as np
import matplotlib.pyplot as plt
from skued import gaussian
image = np.zeros( (256, 256) )
xc, yc = image.shape[0]/2, image.shape[1]/2 # center
extent = np.arange(0, image.shape[0])
xx, yy = np.meshgrid(extent, extent)
rr = np.sqrt((xx - xc)**2 + (yy-yc)**2)
image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
image[np.logical_and(rr < 40, rr > 38)] = 1
image[np.logical_and(rr < 100, rr > 98)] = 0.5
image /= image.max()        # Normalize max to 1
image += np.random.random(size = image.shape)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image)

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()
```

就可以很简单的计算一个角度平均

```py
from skued import azimuthal_average
radius, intensity = azimuthal_average(image, (xc, yc))
plt.plot(radius, intensity)
```

### 4.5 Bonus:去除热点

基线移除（在基线教程中描述）的一个有趣用例是从图像中移除热点。

考虑以下衍射图案：

![image-8](http://scikit-ued.readthedocs.io/en/master/_images/image-8.png)

```py
import matplotlib.pyplot as plt
from skued import diffread

im = diffread('hotspots.tif')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(im, vmin = 0, vmax = 2e3)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()
```

去除热点的代码：

```py
from skued import diffread, baseline_dwt

im = diffread('hotspots.tif')
denoised = baseline_dwt(im, max_iter = 250, level = 1, wavelet = 'sym2', axis = (0, 1))
```

尝试不同组合的小波、级别和迭代次数(max_iter)。

请注意，在图像的情况下使用的基线去除函数是`baseline_dwt()`，它适用于二维数组。 `baseline_dt()`是不可能的，它只适用于一维数组。


## 5 Plotting Tutorial

时间分辨衍射测量可以从精准绘图，特别是时间连续的数据中受益。 本教程将介绍一些`skued`里使其更容易使用的功能。

### color

**time-order**

```py
from skued import spectrum_colors
colors = spectrum_colors(5)
# Colors always go from purple to red (increasing wavelength)
purple = next(colors)
red = list(colors)[-1]
```






# Reference/API

## 1 Baseline-determination

### 1.1 skued.baseline_dt

```py
skued.baseline_dt(array, max_iter, level=None, first_stage='sym6', wavelet='qshift1', background_regions=[], mask=None, mode='constant', axis=-1)
```

基于双树复小波变换的基线确定迭代方法。 

该功能仅适用于1D。

参数：
        - array (~numpy.ndarray, shape (M,N))  背景数据，可以是一维数据或者二维矩阵
        - max_iter (int) 迭代次数