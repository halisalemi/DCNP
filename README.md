# DCNP
This code accompanies the paper "Solving the distance-based critical node problem" and is written in C++. If you want to use or cite this code, please cite the paper.


Compiling the code
--------------------
#### How to run the code for "unweighted instances"

1.0. Download all .cpp and .h files. 

2.0. Start a new C++ project and name it "DCNP". Add all .cpp files to "DCNP" Source files and all .h files to its Header.   

3.0. Build the solution (CTRL+SHIFT+B). 

4.0. Open Windows Command Prompt.

5.0. Find "Release" folder in the "DCNP" folder. Type " cd [folder address of "Release"] " and press Enter. 

6.1. To see the results for "Recursive (R) formulation", type " DCNP.exe Recursive dimacs [folder address of the instance] \ [name of instance].graph [value of k] [value of b] " and press Enter.

6.2. To see the results for "Path-like formulation (PATH)" where k=3, type " DCNP.exe Path_like_k3 dimacs [folder address of the instance] \ [name of instance].graph [value of b] " and press Enter.

6.3. To see the results for "Path-like formulation (PATH)" where k=4, type " DCNP.exe Path_like_k4 dimacs [folder address of the instance] \ [name of instance].graph [value of b] " and press Enter.

6.4. To see the results for "Thin formulation (THIN_I)" using integer separation, type " DCNP.exe Thin_I dimacs [folder address of the instance] \ [name of instance].graph [value of k] [value of b] " and press Enter.

6.5. To see the results for "Thin formulation (THIN_F)" using fractional separation, type " DCNP.exe Thin_F dimacs [folder address of the instance] \ [name of instance].graph [value of k] [value of b] " and press Enter.


#### How to run the code for "weighted instances"

1.0. Download all .cpp and .h files. 

2.0. Start a new C++ project and name it "DCNP". Add all .cpp files to "DCNP" Source files and all .h files to its Header.   

3.0. Build the solution (CTRL+SHIFT+B). 

4.0. Open Windows Command Prompt.

5.0. Find "Release" folder in the "DCNP" folder. Type " cd [folder address of "Release"] " and press Enter.

6.0. Type " DCNP.exe Thin_weighted weighted_graph [folder address of the instance] \ [name of instance].txt [value of k] [value of b] " and press Enter. 


Terms and conditions
--------------------
Copyright (c) 2020 Hosseinali Salemi, Austin L. Buchanan. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
