# Genetic Algorithms Homework

# Task

## Homework T0
Implementati doua metode (una euristica, cealalta determinista) pentru gasirea punctului de maxim sau de minim al unei functii cu un numar arbitrar de variabile.
Testati metoda pe una sau mai multe din functiile de aici.
Intocmiti un raport (format text sau, preferabil, LaTeX + .pdf) in care sa precizati, pentru fiecare functie folosita pentru testare: timpul minim, mediu si maxim de executie cea mai buna si cea mai slaba solutie, precum si media solutiilor obtinute dupa un numar de rulari. Este esential sa faceti o comparatie intre metode, si sa explicati diferentele si rezultatele obtinute. Pana in laboratorul 2. Punctajul, impartit la 10, va fi adaugat punctajului de laborator, ca bonus.

## Homework T1
Find the minimum for the following functions: De Jong 1, Schwefel's, Rastrigin's, Michalewicz's
using the Hill Climbing (both the first improvement and best improvement variants) and Simulated Annealing algorithms.
Soft deadline: lab 5 (-10% points for each week of delay)
Send the homework: source code without binary object and .txt or .pdf file, in a .zip archive (no other archive format allowed) to my e-mail, with a subject written like: [GA]_Name_Surname_Group_T1

## Homework T1'
For the function f=x3-60x2+900x+100, with x in [0, 31], find the maximum.
On [0, 31], the function is uni-modal, having a single maximum point, for which the value is 10.
Study the maximisation of the function using a Hill Climbing algorithm where a candidate solution is represented on 5 bits (32 possible values, so all the integers from 0 to 31). A candidate's neighbourhood is all the bitstrings at a Hamming distance of 1.
Study and explain the function landscape in the context of the first improvement and best improvement variants of the Hill Climbing algorithm. Specify the attraction basing of all local maximum points (attraction basin: the set of points for which the gradient search leads to the same optimum).
You can solve this on paper, +photograph, check the quality of the photo, then send it on e-mail.
Hard deadline: lab 6
Homework T2: Solve Homework T1 using a Genetic Algorithm.
Soft deadline: lab 9

## Homework T2
Solve Homework T1 using a Genetic Algorithm.

## Homework T2'
Optimise the Genetic Algorithm used for Homework T2 by adding adaptive parameters, adaptive operators, Gray coding, etc. Write a short report, detailing and explaining the performance differences.
Hard deadline: lab 10

## Homework T3
Solve one of the problems below with a Genetic Algorithm and another heuristic:
1. The Quadratic Assignment Problem.
2. The Graph Colouring Problem.(Test instances)
3. The Traveling Salesperson Problem (TSP)(Test instances)
4. The Logical Satisfiability Problem (SAT). (Test instances)
Deadline: Second-to-last lab, try to finish it earlier
