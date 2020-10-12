def solver_hook(solve_expr, solve_var, solver_id, guess):
    from dagrt.expression import match, substitute

    pieces = match("unk + (-1)*<func>impl_y(y=sub_y + coeff*unk, t=t)", solve_expr,
                   pre_match={"unk": solve_var})
    pieces["guess"] = guess
    return substitute("<func>solver(sub_y, coeff)", pieces)


def am_solver_hook(solve_expr, solve_var, solver_id, guess):
    from dagrt.expression import match, substitute

    pieces = match(
            "unk + (-1)*<func>f(t=t, fast=sub_fast + coeff*unk, "
            "slow=sub_slow)", solve_expr, pre_match={"unk": solve_var})
    pieces["guess"] = guess
    return substitute("<func>solver(sub_fast, sub_slow, "
                      "coeff, t)", pieces)
