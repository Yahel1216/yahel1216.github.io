---
title: "Clarkson's Algorithm for Linear Programming"
date: 2026-03-31
slug: clarksons-lp-algorithm
draft: false
katex: true
description: "Clarkson's randomized algorithm reduces solving an n×d LP to O(d log(n/d)) sub-problems of size O(d²)×d by reweighting violated constraints — a multiplicative-weights argument in the geometry of basic feasible solutions."
tags: ["linear-programming", "randomized-algorithms", "optimization", "combinatorics"]
categories: ["Optimization", "Algorithms"]
---

# Clarkson's Algorithm for Linear Programming

> [!info] Prerequisites
> *   **Linear Programming:** Feasibility, optimality, and the standard form of an LP.
> *   **Linear Algebra:** Rank, linear independence, matrix inverses, and spanning sets.
> *   **Probability:** Expectation and Markov's inequality.

---

Linear programs arise throughout combinatorial optimization, machine learning, and operations research. In many practical settings the number of constraints $n$ is enormous compared to the ambient dimension $d$ — think of $n = 10^6$ constraints in $d = 50$ dimensions. Standard algorithms like the simplex method take time proportional to $n$ per pivot, and interior-point methods scale as $O(n^{3.5})$ in the worst case. When $n \gg d$, most constraints are redundant: the optimal solution is determined by only $d$ of them.

In this post we cover **Clarkson's randomized algorithm** [[1]](#references), which exploits this redundancy. The algorithm reduces the original $n \times d$ LP to $O(d \log(n/d))$ sub-problems, each involving only $O(d^2)$ constraints. The key engine is a *multiplicative-weights* reweighting scheme: constraints that are frequently violated accumulate higher weight, rapidly concentrating the important structure. In this post, we develop all the geometric prerequisites — extreme points, basic feasible solutions, and the connection between the two — and then prove the main runtime bound from scratch.

## Some Intuition

Recall that a linear program is an optimization problem of the form:

$$\left\lbrace \begin{array}{ll} \mathrm{minimize} & c^{\top}x \\\\ \mathrm{subject\ to} & Ax\le b \\\\ & x\in\mathbb{R}^{d} \end{array}\right.$$

where $A\in\mathbb{R}^{n\times d}$ is the constraints matrix, $b\in\mathbb{R}^{n}$ the bias, and $c\in\mathbb{R}^{d}$ the objective. The **feasible set**

$$\mathcal{P}=\lbrace x\in\mathbb{R}^{d}:Ax\le b\rbrace$$

is a convex set called a **polytope**. We assume it is bounded: there exists $R>0$ such that $\mathcal{P}\subset\lbrace x\in\mathbb{R}^{d}:\|x\|\le R\rbrace$.

### Extreme Points

Intuitively, a polytope is a convex body with "corners." These corners are the extreme points:

> [!caution] Definition
> For a convex set $C$, a point $x\in C$ is an **extreme point** if:
> $$\forall\, y,z\in C,\ t\in (0,1):\quad x=t\cdot y+(1-t)\cdot z\implies x=y=z.$$
> In other words, a point is extreme if it cannot be written as the (non-trivial) convex combination of two distinct points.

For a point $x\in\mathcal{P}$, define its **active set** $B(x)\subset [n]$ as the constraints that are tight at $x$:

$$B(x) = \lbrace i \in [n] : a_i^{\top}x = b_i\rbrace,$$

where $a_i$ is the $i$-th row of $A$. The active set captures exactly which directions $x$ is limited at, meaning that every trajectory $v$ that satisfies $a_i^{\top}v>0$ for some $i\in B(x)$ is not a valid trajectory to move from $x$. That is, $x+\epsilon v$ is infeasible, no matter how small $\epsilon>0$ we choose.

> [!important] Lemma
> A point $x\in \mathcal{P}$ is an extreme point if and only if the active set $B(x)$ satisfies $|B(x)|\ge d$ and the set $\lbrace a_{i}:i\in B(x)\rbrace$ spans $\mathbb{R}^{d}$.

**Proof.** Denote $B=B(x)$ and $U=\mathrm{Span}\lbrace a_{i}:i\in B\rbrace$. Note that $\lbrace a_i : i \in B \rbrace$ spanning $\mathbb{R}^d$ automatically implies $|B| \geq d$, so we only need to check the spanning condition.

Suppose $x$ is an extreme point and assume for contradiction that $\dim U < d$. In this case $\dim U^{\perp}\ge 1$ (the orthogonal complement of $U$). Let $0\ne y\in U^{\perp}$ and set $z=x+\varepsilon y$ for $\varepsilon>0$ to be chosen later. For $i\in B$, since $y \perp a_i$, we have $$a_{i}^{\top}z=a_{i}^{\top}x+\varepsilon\cdot a_{i}^{\top}y=b_{i}.$$
For $i\notin B$, we have $a_{i}^{\top}x < b_{i}$ by definition, so
$$\varepsilon_{0}=\min_{i\notin B}(b_{i}-a_{i}^{\top}x)>0,$$
since this is a finite minimum. Setting $M=\max_{i\notin B}|a_{i}^{\top}y|$, if $\varepsilon< \varepsilon_{0}/M$ then for $i\notin B$:

$$b_{i}-a_{i}^{\top}z=b_{i}-a_{i}^{\top}x-\varepsilon\cdot a_{i}^{\top}y\ge \varepsilon_{0}-\varepsilon M>0,$$

so $z$ is feasible. The same argument shows $w=x-\varepsilon y$ is feasible. But then $x=\frac{1}{2}z+\frac{1}{2}w$ with $z\ne w$ (since $y\ne 0$), contradicting that $x$ is extreme.

Conversely, suppose $B$ spans $\mathbb{R}^{d}$. The system $A_{B}x=b_{B}$ (restricting rows to $B$) has $x$ as a solution, and since $A_B$ has full rank the solution is **unique**. Suppose $x=ty+(1-t)z$ for $y,z\in\mathcal{P}$ and $t\in(0,1)$. Then $A_{B}y\le b_{B}$ and $A_{B}z\le b_{B}$, and
$$b_{B}=A_{B}x=t\cdot A_{B}y+(1-t)\cdot A_{B}z\le t\cdot b_{B}+(1-t)\cdot b_{B}=b_{B}.$$
Equality forces $A_{B}y=A_{B}z=b_{B}$ and by uniqueness, we conclude that $x=y=z$. Therefore $x$ is an extreme point by definition. $\blacksquare$

### From Extreme Points to Vertices

A common and convenient assumption is that the active set is exactly of size $d$ at each extreme point.

> [!caution] Definition
> In a polytope $\mathcal{P}=\lbrace x:Ax\le b\rbrace$, extreme points are called **vertices**. The polytope is **non-degenerate** if for every vertex $x$, the active set $B(x)$ has exactly size $d$.

The standard technique to pass to the non-degenerate case is a small random perturbation: for a Gaussian random vector $\xi\in\mathbb{R}^{n}$, the perturbed polytope $\mathcal{P}'=\lbrace x:Ax\le b+\xi\rbrace$ is non-degenerate with probability 1. For small enough $\xi$, the two polytopes define essentially the same problem. **From now on we assume $\mathcal{P}$ is non-degenerate.**

> [!caution] Definition
> Let $B\subset [n]$ with $|B|=d$ such that $A_{B}$ is invertible. The **basic solution** associated with $B$ is $x=A_{B}^{-1}b_{B}$. If this solution also satisfies $Ax\le b$, it is a **Basic Feasible Solution (BFS)** and $B$ is called a **basis** for $x$.

> [!important] Lemma
> A point $x\in \mathcal{P}$ is a vertex if and only if it is a BFS.

**Proof.** By non-degeneracy and the preceding lemma, $x$ is a vertex iff $B(x)$ has size $d$ and spans $\mathbb{R}^{d}$, i.e., iff $A_{B(x)}$ is invertible. Then $x=A_{B(x)}^{-1}b_{B(x)}$ is a BFS with basis $B(x)$. Conversely, if $x$ is a BFS with basis $B$, then $B\subseteq B(x)$, so $B(x)$ spans $\mathbb{R}^{d}$, and $x$ is a vertex. $\blacksquare$

Vertices correspond bijectively to feasible bases. There are at most $\binom{n}{d}$ bases, and the feasible ones form a subset. Since $n \gg d$, this is a vast search space.

### The Optimum is Always a Vertex

The Krein–Milman theorem states that any compact convex set is the convex hull of its extreme points. Denoting by $\mathcal{B}(A)$ the set of all BFSs:

> [!important] Corollary
> Every point in $\mathcal{P}$ is a convex combination of BFSs:
> $$\mathcal{P}=\left\lbrace \sum_{j=1}^{m}\alpha_{j}x_{j} \middle| \sum_{j=1}^{m}\alpha_{j}=1,\ \alpha_{j}\ge 0,\ x_{j}\in \mathcal{B}(A)\right\rbrace.$$

> [!important] Lemma
> The optimum of the LP is attained at a BFS.

**Proof.** For any $y\in\mathcal{P}$, write $y=\sum_j\alpha_j x_j$ with $x_j\in\mathcal{B}(A)$, $\alpha_j\ge 0$, $\sum_j\alpha_j=1$. Then $$c^{\top}y=\sum_j\alpha_j\cdot  c^{\top}x_j\le\max_j c^{\top}x_j$$using the non-negativity of the coefficients and the fact they sum to $1$. Hence $\max_{y\in\mathcal{P}}c^{\top}y=\max_{x\in\mathcal{B}(A)}c^{\top}x.$ $\blacksquare$

By a small perturbation of $c$ we can also assume the optimal BFS is **unique**. If there is a tie between two BFSs with different bases, a small perturbation of $c$ will break the tie.

The goal of LP is therefore to identify the basis of the optimal BFS. With $n \gg d$, the search space is huge. The key insight of Clarkson's algorithm is that we can rapidly identify which constraints are important by reweighting them according to how often they are violated.

## The Algorithm

We work with **multi-sets**: collections in which repetitions are allowed. For a multi-set $S$, write $|S|$ for its cardinality (counting multiplicity), and $\binom{S}{r}$ for the collection of all sub-multi-sets of size $r$. For example the cardinality of the multi-set $\lbrace 1,1,2,3\rbrace$ is $4$.

The heart of the algorithm is the following probabilistic estimate.
> [!important] Lemma (Low Violation Lemma)
> Let $S$ be a multi-set of constraints, $r \geq d$, and $R\sim\binom{S}{r}$ chosen uniformly at random. Let $x_{R}^{\*}$ be the optimal solution to
> $$\mathbf{P}\_{R}:\quad \mathrm{minimize}\ c^{\top}x\quad \mathrm{subject\ to}\ a_{i}^{\top}x\le b_{i}\ \forall i\in R,$$
> and let $V_{R}=\lbrace i\in S: a_{i}^{\top}x_{R}^{\*}>b_{i}\rbrace$ be the constraints of $S$ violated by $x_R^{\*}$. Then:
> $$\mathbb{E}\_{R\sim\binom{S}{r}}[|V_{R}|]\le d\cdot \frac{|S|-r}{r+1}.$$

**Proof.** For $R\in\binom{S}{r}$ and $i\in S$, let $v(i,R)=1$ if $i\in V_R$ and $0$ otherwise. Note $v(i,R)=0$ whenever $i\in R$, since $x_R^{\*}$ is feasible for $\mathbf{P}\_R$. We can write
$$\binom{|S|}{r}\cdot \mathbb{E}\_{R}[|V_R|]=\sum\_{R\in\binom{S}{r}}|V_{R}|=\sum\_{R\in\binom{S}{r}}\sum_{i\in S\setminus R}v(i,R)=\sum_{Q\in\binom{S}{r+1}}\sum_{i\in Q}v(i,Q\setminus\lbrace i\rbrace),$$

since choosing $R$ and then $i\in S\setminus R$ is equivalent to choosing $Q$ of size $r+1$ and removing one element $i\in Q$.

Fix $Q\in\binom{S}{r+1}$, and let $B$ be the basis of the unique optimal solution $x_{Q}^{\*}$ to $\mathbf{P}\_{Q}$. For $i\in Q$ to have $v(i,Q\setminus\lbrace i\rbrace)=1$, meaning $x_{Q\setminus\lbrace i\rbrace}^{\*}$ violates $i$, it is necessary that $i\in B$. If instead $i\notin B$, then $B\subset Q\setminus\lbrace i\rbrace$, so by uniqueness $x_{Q\setminus\lbrace i\rbrace}^{\*}=x_{Q}^{\*}$, and $x_{Q}^{\*}$ satisfies all constraints in $Q$ including $i$. Therefore $\sum_{i\in Q}v(i,Q\setminus\lbrace i\rbrace)\le|B|=d$, and:

$$\mathbb{E}[|V_{R}|]=\frac{\sum_{R}|V_R|}{\binom{|S|}{r}}\le\frac{d\cdot\binom{|S|}{r+1}}{\binom{|S|}{r}}=d\cdot\frac{|S|-r}{r+1}.\quad\blacksquare$$

The bound says that no matter how $S$ is changed, the expected number of violations is at most a constant fraction $d/(r+1)$ of $S$. This is the invariant the algorithm exploits.

> [!caution] Algorithm: Clarkson's LP Algorithm
> *   **Input:** A matrix $A\in\mathbb{R}^{n\times d}$, a vector $b\in\mathbb{R}^{n}$, and an objective $c\in\mathbb{R}^{d}$.
> *   **Output:** The optimal BFS for the LP.
>
> 1.  Initialize $S\gets [n]$ as a multi-set (each constraint with multiplicity 1).
> 2.  Set $r = C d^{2} \log e$ for some constant $C > 2$ (see the analysis).
> 3.  Repeat:
>     1.  Draw $R\sim\binom{S}{r}$ uniformly at random.
>     2.  Solve $\mathbf{P}\_{R}$ to obtain the optimal BFS $x_{R}^{\*}$.
>     3.  Compute $V_{R}=\lbrace i\in S: a_{i}^{\top}x_{R}^{\*}>b_{i}\rbrace$.
>     4.  If $V_{R}=\emptyset$, return $x_{R}^{\*}$.
>     5.  If $|V_{R}|\le 2d\cdot\dfrac{|S|-r}{r+1}$: append $V_{R}$ to $S$, meaning $S\gets S+V_{R}$. This multiplies the violations' multiplicities in $S$.

The condition in step 3.5 checks whether $|V_R|$ is at most twice its expected value (the factor of 2 comes from Markov's inequality, as we will see). When satisfied, the violated constraints are reweighted by doubling their multiplicity in $S$.

## The Runtime

The algorithm terminates when $V_R = \emptyset$ (step 3.4), which means $x_R^{\*}$ is feasible for every constraint in $S$ — and since $S$ always contains all original constraints, $x_R^{\*}$ is a globally optimal BFS.

By Markov's inequality:

$$\Pr \left[|V_{R}|> 2\cdot\mathbb{E}[|V_{R}|]\right]\le\frac{1}{2},$$

so the condition in step 3.5 is satisfied with probability at least $\frac{1}{2}$. The expected number of total loop iterations is at most twice the number of **successful** iterations (those satisfying the condition). We now bound the latter.

> [!important] Lemma
> The condition $|V_{R}|\le \rho|S|$ can be satisfied $kd$ times, as long as $k$ satisfies:
> $$k(1-d\rho\log e)\le \log n-\log d.$$
> Here $\log$ denotes $\log_2$.

**Proof.** Index only the successful iterations as $t=1,2,\ldots$, and let $S_t$ denote the state of $S$ at the end of iteration $t$. Let $B^{\*}$ be the basis for the globally optimal solution. For each $j\in B^{\*}$, write $m_j^t$ for the multiplicity of $j$ in $S_t$, and let $B_t$ denote the corresponding sub-multi-set of $S_t$.

Initially $m_j^0=1$ for all $j$, and $|B_0|=d$. Note that if $V_R\not=\emptyset$, then $V_R \cap B^{\*}$ is not empty. Otherwise, $x_R^\*$ is also an optimal solution to the problem $\mathbf{P}\_{R\cup B^\*}$, whose feasible set is contained in that of $\mathbf{P}\_{R}$. Since the optimal BFS $x^\*$ is the unique optimal solution to $\mathbf{P}\_{S}\subset \mathbf{P}\_{R\cup B^*}\subset \mathbf{P}\_{B^\*}$, it follows that $x^\*=x_R^\*$ (in this notation we mean the feasible sets are contained).

Now, note that if $j\in B^\*$ is violated, then all $m_j^t$ copies of $j$ in $S_t$ appear in $V_R$. Thus adding $V_R$ to $S$ doubles the multiplicity of at least one basis constraint. After $t=id$ iterations, the cardinality of $|B_t|$ must be at least $d\cdot 2^i$. On the other hand, since $|V_R|\le \rho |S_t|$, the cardinality of $S$ grows by a factor of $(1+\rho)$ at most. After $id$ iterations, it holds $|S_{id}|\le (1+\rho)^{id}\cdot |S_0|\le e^{id\rho} \cdot n$. Since $B_{id}\subset S_{id}$, we obtain $d\cdot 2^i\le e^{id\rho}\cdot n$. Taking logarithms,
$$\log d+i\le id\rho\log e+\log n,$$
which rearranges to $i(1-d\rho\log e)\le\log(n/d)$. $\blacksquare$

### Choosing $r$

The algorithm's success condition $|V_R|\le 2d\cdot\frac{|S|-r}{r+1}$ is exactly the requirement $|V_R|\le\rho|S|$ with $\rho=\frac{2d(|S|-r)}{(r+1)|S|}\le\frac{2d}{r+1}$, so we take $\rho=\frac{2d}{r+1}$. The choice of $r$ creates a fundamental trade-off:

- **Smaller $r$:** Sub-problems $\mathbf{P}_R$ are cheaper to solve, but violations are more numerous, requiring more reweighting steps.
- **Larger $r$:** Fewer reweighting steps, but each sub-problem involves more constraints.

Three cases illustrate the trade-off:

1. **$r=d^{2}$.** Then $d\rho\log e\approx 2\log e>1$, so $1-d\rho\log e<0$. The bound gives no useful upper bound on $k$, and the algorithm is not guaranteed to terminate.

2. **$r=Cd^{2}\log e$ for $C>2$.** Then $\rho\le\frac{2}{Cd\log e}$, so $d\rho\log e=\frac{2}{C}<1$. The bound becomes $k\left(1-\frac{2}{C}\right)\le\log(n/d)$, giving $k\le\frac{C\log(n/d)}{C-2}$. The expected total successful iterations is $kd=O(d\log(n/d))$.

3. **$r=2d^{2+\alpha}\log e$ for $\alpha>0$.** Then $d\rho\log e=\frac{1}{d^{\alpha}}$, giving $k\le\frac{d^{\alpha}\log(n/d)}{d^{\alpha}-1}$. The iteration count improves by only a constant factor over case 2, while each sub-problem is much more expensive.

The sweet spot is $r=\Theta(d^{2})$: **the algorithm terminates after $O(d\log(n/d))$ successful reweighting steps on average, and each step solves a sub-problem with $O(d^{2})$ constraints in $d$ dimensions.**

### Summary of the Complexity

We have reduced a single $n\times d$ LP to $O(d\log(n/d))$ sub-problems of size $O(d^{2})\times d$. When $n\gg d$, this is a dramatic improvement: for instance, $n=10^6$ constraints in $d=100$ dimensions gives roughly a few hundred sub-problems, each involving only $10^4$ constraints.

### Where Is the Analysis Loose?

The tradeoff may not be optimal. Two sources of slack are worth noting.

First, the lower bound on $B_{t}$ assumes that at each iteration only a single basis element doubles. In practice, $V_R$ may contain multiple basis elements, causing much faster growth. This may be affected by properties like conditioning of the constraints matrix (measuring similarity between rows).

Second, the bound $|S_t|\le|S_{t-1}|(1+\rho)$ assumes $|V_R|$ is always near its maximum of $\rho|S|$. If $|V_R|$ is typically much smaller — as one would expect when the constraints are structured — the growth of $|S|$ slows considerably and the algorithm converges faster.

Let's see how progress on one of these may affect the runtime:
> [!important] Lemma
> Let $a_t$ denote the number of distinct basis constraints violated at $t$-th (successful) iteration. Suppose the algorithm runs for $T$ successful iterations, then $T\le O\left(\frac{d}{a} \log (n/d)\right)$, where $a=\frac{1}{T}\sum_{t}a_t$ is the average.

**Proof.** This is a very similar proof as before, but formalized with a potential function. Define $\phi_t= \prod_{j\in B^{\*}}m_j^t$ (the product of multiplicities of basis elements). Then $\phi_T \ge 2^{\sum_t a_t}=2^{Ta}$, as this is the minimal way to double the required number of elements. By the arithmetic-geometric means inequality, $\sqrt[d]{\phi_T} \le \frac{\sum_{j\in B^\*}m_j^T}{d}= \frac{|B_{T}|}{d}$, which gives $2^{Ta/d}\le |S_T|/ d\le e^{T\rho}n/d$. Therefore $d 2^{Ta/d} \le e^{T\rho}n$ which implies $\frac{Ta}{d}-T\rho \log e\le \log n -\log d$, that is $T(a/d-\rho\log e)\le \log(n/d)$. Substituting $\rho\le \frac{a}{2d\log e}$ we obtain $T\le O\left( \frac{d}{a}\log (n/d)\right)$. This setting of $\rho$ is achieved for $r\ge 4\log e\cdot d^2 / a$, which can be much smaller than $d^2$. $\blacksquare$

Whether $r=\Theta(d^2)$ can be improved remains an interesting open question in the general case. Tighter results may be possible under structural assumptions on $A$. A possible way to prove such a result is a **win-win** argument. For example, proving that the average $a$ (as defined above) is rather large, OR some other useful property must hold.

## Conclusion

Clarkson's algorithm is a beautiful instance of the **multiplicative weights method** applied to a geometric problem. The key facts that make it work are: (1) the optimal solution is always a vertex, determined by only $d$ of the $n$ constraints; (2) a random sample of $O(d^2)$ constraints already "sees" the optimal basis with reasonable probability, as quantified by the Low Violation Lemma; and (3) reweighting violated constraints rapidly concentrates the sampling distribution near the important constraints.

The result is a reduction from an $n \times d$ LP to $O(d \log(n/d))$ problems of size $O(d^2) \times d$ — an algorithm that is both practically fast and theoretically clean, almost independent of $n$, which exploits the geometry of linear programming in a fundamental way.

## References

1.  **Clarkson, K.L.** (1995). [*Las Vegas Algorithms for Linear and Integer Programming When the Dimension is Small*](https://doi.org/10.1145/201019.201036).
2.  **Seidel, R.** (1991). [*Small-Dimensional Linear Programming and Convex Hulls Made Easy*](https://link.springer.com/article/10.1007/BF02574699).
3.  **Bertsimas, D., & Tsitsiklis, J.N.** (1997). [*Introduction to Linear Optimization*](https://www.researchgate.net/publication/235558951_Introduction_to_Linear_Optimization). 
4.  **Arora, S., Hazan, E., & Kale, S.** (2012). [*The Multiplicative Weights Update Method: a Meta-Algorithm and Applications*](https://doi.org/10.4086/toc.2012.v008a006).
