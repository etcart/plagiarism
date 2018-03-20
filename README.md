# Plagiarism
Plagiarism detector for student works, based on RIVet system.
The purpose of this system is to be able to identify student submitted C code which is at a high risk of plagiarism.  
This is done by random index vector analysis of underlying assembly code.

1. collect all student works under 1 root directory
2. add to that root directory, online code that code might have been taken from
3. run python script vindicators.py 
4. will print a list of high risk pairs, along with level of similarity  (out of 1)

for details on the underlying algorithm, see https://github.com/etcart/RIVet-C the core library of random index vector tools
