# Adaptive Constraint Screening for MMA-Based Structural Topology Optimization

## Purpose

This document describes a practical approach for reducing expensive constraint sensitivity/adjoint calculations in structural topology optimization with many load cases while keeping the MMA optimizer itself unchanged.

The central idea is:

> Evaluate all required constraint **values**, but compute expensive **sensitivities only for active or near-active constraints**.  
> If a skipped constraint becomes important at the MMA trial design, reject that trial, activate the constraint, compute its sensitivity at the original current design, and rerun MMA from the same design.

---

## 1. Terminology

Let

\[
\mathbf{x}_k
\]

be the **current accepted design** at optimization iteration \(k\).

For density-based topology optimization,

\[
\mathbf{x}_k=[x_1,x_2,\ldots,x_N]^T,
\qquad 0\le x_j\le1.
\]

Let

\[
\mathbf{x}_{trial}
\]

be the **candidate design proposed by MMA**.

It is not yet an accepted optimization iterate.

Only after the trial design passes the required checks do we set

\[
\mathbf{x}_{k+1}=\mathbf{x}_{trial}.
\]

Therefore:

```text
x_k = current accepted design
        |
        v
       MMA
        |
        v
x_trial = proposed design
        |
        v
check all constraints
        |
     +--+--+
     |     |
   reject accept
     |     |
 stay at  x_(k+1) = x_trial
   x_k
```

---

## 2. Main Recommendation

### Do not modify the mathematical MMA algorithm.

Implement a **Constraint Screening Manager outside MMA**:

```text
FEA / structural analysis
        |
        v
Evaluate ALL constraint values
        |
        v
Constraint Screening Manager
        |
        +--> select active / near-active constraints
        |
        v
Compute sensitivities ONLY for selected constraints
        |
        v
Existing MMA solver
        |
        v
x_trial
        |
        v
Evaluate ALL constraints at x_trial
        |
        +--> skipped constraint became important?
                 |
            +----+----+
            |         |
           YES        NO
            |         |
      reject trial    accept trial
      activate it
      rerun MMA
      from same x_k
```

The original constraints remain part of the optimization problem. Screening controls which sensitivities are calculated and which constraint rows are supplied to the current reduced MMA subproblem.

---

## 3. Why Return to the Same x_k?

Suppose the current active set is

```text
C1, C2, C3
```

Their sensitivities are evaluated at the current design:

\[
\nabla C_1(\mathbf{x}_k),\quad
\nabla C_2(\mathbf{x}_k),\quad
\nabla C_3(\mathbf{x}_k).
\]

MMA uses this information to construct its local subproblem around \(\mathbf{x}_k\) and proposes \(\mathbf{x}_{trial}\).

Now suppose the full response check at \(\mathbf{x}_{trial}\) reveals that a previously skipped constraint \(C_{17}\) has become near-active or violated.

This means \(C_{17}\) was relevant to the proposed step.

Do **not** accept the trial.

Instead:

1. Keep the accepted design at \(\mathbf{x}_k\).
2. Add \(C_{17}\) to the active set.
3. Compute

\[
\nabla C_{17}(\mathbf{x}_k).
\]

4. Rerun MMA from the same \(\mathbf{x}_k\).
5. Check the new trial design again.

The sequence is:

```text
x_k
 |
 +--> MMA with C1,C2,C3
 |          |
 |          v
 |      x_trial_1
 |          |
 |     C17 becomes important
 |          |
 |        REJECT
 |
 +--> compute dC17/dx at x_k
 |
 +--> MMA with C1,C2,C3,C17
            |
            v
        x_trial_2
            |
       full check OK
            |
          ACCEPT
            |
            v
      x_(k+1) = x_trial_2
```

The rejected trial is only a proposal. The outer optimization iteration has not advanced.

---

## 4. Recommended Algorithm

```text
for optimization iteration k:

    # --------------------------------------------------
    # 1. Current accepted design
    # --------------------------------------------------

    x = x_k


    # --------------------------------------------------
    # 2. Solve structural problem at x_k
    # --------------------------------------------------

    state_k = solve_FEA(x_k)


    # --------------------------------------------------
    # 3. Evaluate ALL constraint values
    # --------------------------------------------------

    g_all = evaluate_all_constraints(x_k, state_k)


    # --------------------------------------------------
    # 4. Select active / near-active constraints
    # --------------------------------------------------

    active = select_active_constraints(
        g_all,
        previous_active_set
    )


    # --------------------------------------------------
    # 5. Compute objective and objective sensitivity
    # --------------------------------------------------

    f_k    = objective(x_k, state_k)
    df0dx  = objective_sensitivity(x_k, state_k)


    # ==================================================
    # INNER ACTIVE-SET REPAIR LOOP
    # ==================================================

    while True:

        # ----------------------------------------------
        # 6. Compute sensitivities ONLY for active set
        # ----------------------------------------------

        dgdx_active = sensitivity(
            design      = x_k,
            constraints = active,
            state       = state_k
        )


        # ----------------------------------------------
        # 7. Run normal MMA
        # ----------------------------------------------

        x_trial = MMA(
            x           = x_k,
            objective   = f_k,
            df0dx       = df0dx,
            constraints = g_all[active],
            gradients   = dgdx_active,
            state       = original_MMA_state_k
        )


        # ----------------------------------------------
        # 8. Evaluate ALL constraints at trial design
        # ----------------------------------------------

        state_trial = solve_FEA(x_trial)

        g_trial = evaluate_all_constraints(
            x_trial,
            state_trial
        )


        # ----------------------------------------------
        # 9. Find skipped constraints that became
        #    near-active or violated
        # ----------------------------------------------

        newly_active = []

        for constraint i not in active:

            if g_trial[i] >= reactivation_threshold:

                newly_active.add(i)


        # ----------------------------------------------
        # 10. Repair active set if necessary
        # ----------------------------------------------

        if newly_active is not empty:

            active = active UNION newly_active

            # IMPORTANT:
            #
            # Do NOT accept x_trial.
            # Do NOT increment k.
            # Do NOT update accepted MMA history.
            #
            # Newly required sensitivities are computed
            # at x_k on the next pass through this loop.

            continue


        # ----------------------------------------------
        # 11. Candidate passed screening check
        # ----------------------------------------------

        break


    # ==================================================
    # END INNER LOOP
    # ==================================================


    # --------------------------------------------------
    # 12. Accept candidate
    # --------------------------------------------------

    x_(k+1) = x_trial


    # --------------------------------------------------
    # 13. Update MMA history/state ONCE
    # --------------------------------------------------

    update_MMA_state(
        x_k,
        x_(k+1)
    )


    # --------------------------------------------------
    # 14. Near convergence, verify full problem
    # --------------------------------------------------

    if near_convergence:

        disable_screening()

        compute_all_required_current_sensitivities()

        verify_full_problem()


    # --------------------------------------------------
    # 15. Continue or terminate
    # --------------------------------------------------

    if converged:
        return x_(k+1)

    x_k = x_(k+1)
```

---

## 5. Constraint Selection

Normalize constraints whenever practical. For an upper-bound response,

\[
g_i(\mathbf{x})
=
\frac{r_i(\mathbf{x})}
     {r_{i,\mathrm{allow}}}
-1
\le0.
\]

Then:

- \(g_i=0\): exactly at the limit.
- \(g_i>0\): violated.
- \(g_i<0\): satisfied.

A simple initial screening strategy is:

```text
if g_i >= activation_threshold:
    ACTIVE

elif constraint was active recently:
    KEEP ACTIVE

else:
    sensitivity may be skipped
```

### Use Hysteresis

Do not use exactly the same threshold for activation and deactivation.

For example:

```text
Activate:       g_i > -0.10
Deactivate:     g_i < -0.20
```

The region between them prevents repeated switching:

```text
                 KEEP ACTIVE
                      |
                      |
--------+-------------+-------------+---->
      -20%          -10%            0
        |                           |
   may deactivate              constraint limit

        <---- hysteresis ---->
```

The 10% and 20% values are examples, not universal constants. They should be tuned to the problem, normalization, MMA move limits, and observed constraint variation.

---

## 6. What MMA Receives

Suppose the complete optimization problem contains 60 constraints:

```text
50 displacement constraints
10 stress constraints
```

At iteration \(k\), suppose screening selects

```text
active = [C2, C5, C11, C27, C31]
```

Then MMA receives only:

```text
m = 5

constraint values:
    C2
    C5
    C11
    C27
    C31

constraint gradients:
    dC2/dx
    dC5/dx
    dC11/dx
    dC27/dx
    dC31/dx
```

The other 55 constraints are **not supplied as fake rows** to that reduced MMA call.

They still exist in the global problem registry and are checked at the trial design.

---

## 7. What NOT to Do

### Do not zero skipped gradients

Do not do this:

```text
if inactive:
    dC_i/dx = 0
```

A zero derivative is generally not the true derivative and can create an incorrect MMA approximation.

Instead, omit that constraint from the current reduced MMA subproblem.

---

### Do not accept x_trial before checking skipped constraints

Wrong:

```text
x_k
 -> MMA
 -> x_trial
 -> immediately set x_(k+1) = x_trial
```

Correct:

```text
x_k
 -> MMA
 -> x_trial
 -> evaluate all constraints
 -> accept/reject
```

---

### Do not compute a newly activated gradient at x_trial when repairing the same MMA step

If \(C_{17}\) becomes important at the trial point, calculate

\[
\nabla C_{17}(\mathbf{x}_k),
\]

not

\[
\nabla C_{17}(\mathbf{x}_{trial}),
\]

when repairing the subproblem based at \(\mathbf{x}_k\).

---

### Do not advance the optimization iteration after a rejected trial

During active-set repair:

```text
k stays unchanged
x_k stays unchanged
accepted MMA history stays unchanged
```

Only after acceptance:

```text
x_(k+1) = x_trial
k -> k+1
```

---

## 8. Where the Computational Saving Comes From

Without screening:

```text
x_k
 |
 +--> solve FEA
 |
 +--> evaluate all constraint values
 |
 +--> compute ALL constraint sensitivities/adjoints
 |        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 |                   expensive
 |
 +--> MMA
```

With screening:

```text
x_k
 |
 +--> solve FEA
 |
 +--> evaluate all constraint values
 |
 +--> identify important constraints
 |
 +--> compute sensitivities/adjoints ONLY for them
 |        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 |                   main saving
 |
 +--> MMA
 |
 +--> full response check
```

The primary saving is therefore in **adjoint solves and sensitivity assembly**.

The method may require an additional forward analysis at a trial design. It is most useful when sensitivity/adjoint calculations are expensive and the screening rule predicts important constraints well enough that trial rejection is not excessive.

---

## 9. Additional Exact Optimization Before Screening

Before relying heavily on constraint screening, exploit exact sensitivity reuse wherever possible.

For linear elasticity,

\[
K(\mathbf{x})u_l=f_l.
\]

If many displacement constraints across different load cases use the same displacement measurement operator, their adjoint equation can share the same right-hand side/operator structure.

Therefore, multiple constraints may require fewer independent adjoint solves than the raw number of constraints suggests.

Also:

- reuse the stiffness matrix factorization across load cases when applicable;
- use multiple-right-hand-side solves;
- exploit self-adjoint compliance sensitivities;
- identify identical or linearly dependent adjoint right-hand sides.

These improvements are attractive because they reduce computational cost **without screening out constraints**.

---

## 10. Recommended Software Architecture

Keep the components separated:

```text
+-----------------------------+
| FEA / Structural Solver     |
+-------------+---------------+
              |
              v
+-----------------------------+
| Response Evaluator          |
| - evaluate ALL constraints  |
+-------------+---------------+
              |
              v
+-----------------------------+
| Constraint Screening Manager|
| - active set                |
| - hysteresis                |
| - reactivation              |
| - permanent constraint IDs  |
+-------------+---------------+
              |
              v
+-----------------------------+
| Adjoint/Sensitivity Solver  |
| - selected constraints only |
+-------------+---------------+
              |
              v
+-----------------------------+
| Existing MMA                |
| - unchanged mathematics     |
+-------------+---------------+
              |
              v
          x_trial
              |
              v
+-----------------------------+
| Full Response Check         |
+-------------+---------------+
              |
       +------+------+
       |             |
    reject         accept
       |             |
       v             v
 repair set       x_(k+1)
```

This keeps MMA generic and puts adaptive constraint management in the optimization driver.

---

## 11. Final Recommended Strategy

Use the following hierarchy:

### First: exact cost reduction

Exploit shared adjoints, shared matrix factorizations, multiple RHS solves, and self-adjoint responses.

### Second: external constraint screening

Evaluate all constraint values but calculate expensive sensitivities only for active/near-active constraints.

### Third: trial-design safety check

After MMA proposes \(\mathbf{x}_{trial}\), evaluate all original constraints.

### Fourth: active-set repair

If a skipped constraint becomes important:

```text
reject x_trial
        |
        v
keep x_k
        |
        v
activate constraint
        |
        v
compute its sensitivity at x_k
        |
        v
rerun MMA from x_k
```

### Fifth: full verification near convergence

Disable screening near convergence and verify the complete original optimization problem using all required current constraint information.

---

## 12. Important Mathematical Qualification

This workflow is a **practical screened-MMA method**.

It keeps all original constraints as design requirements and provides a mechanism for reactivating constraints that become important. However, simple threshold-based screening does **not guarantee the identical iteration sequence or identical MMA subproblem solution** that would be obtained by supplying every constraint gradient at every iteration.

If exact equivalence to full-constraint MMA is required, omitted constraints need a mathematically certified redundancy/screening test for the MMA subproblem, or their derivatives must be obtained through exact reuse rather than omitted.

For most implementations whose objective is to reduce expensive adjoint work while retaining the original constraints and the existing MMA solver, the external screening-and-reactivation architecture above is the recommended starting point.