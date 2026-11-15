import Lean

/-!
A deliberately restricted, source-reading semantics for the native primitive.
Only the eight-statement shape below is accepted. Arithmetic expressions are
parsed, not supplied as a hand-written replacement for the source expression.
Python integer // and % are supported only with a positive divisor. List copy,
ordered augmented updates and the returned coin/current are explicit.
This is not a semantics or verified compiler for general Python.
-/

namespace EventAlgebra.NativeBinarySource

inductive Expr where
  | lit : Int → Expr
  | var : String → Expr
  | loadA | loadB
  | abs : Expr → Expr
  | bin : String → Expr → Expr → Expr
  deriving Repr, DecidableEq

structure Program where
  difference : Expr
  current : Expr
  subtractA : Bool
  updateA : Expr
  subtractB : Bool
  updateB : Expr
  recordCoin : Expr
  recordCurrent : Expr
  deriving Repr, DecidableEq

def compact (s : String) : String :=
  String.ofList (s.toList.filter fun c => !c.isWhitespace)

-- These helpers are structural, so the kernel can reduce the exact source.
-- The general String splitOn/toInt?/startsWith/endsWith APIs are intentionally
-- absent from the proof path.
def splitChars (separator : Char) : List Char → List (List Char)
  | [] => [[]]
  | c::cs =>
    let rest := splitChars separator cs
    if c == separator then []::rest
    else match rest with
      | [] => [[c]]
      | first::others => (c::first)::others

def splitText (separator : Char) (s : String) : List String :=
  (splitChars separator s.toList).map String.ofList

def starts (stem s : String) : Bool :=
  s.toList.take stem.length == stem.toList

def ends (stem s : String) : Bool :=
  s.toList.reverse.take stem.length == stem.toList.reverse

def decimal (s : String) : Option Int :=
  let chars := s.toList
  if chars.isEmpty || (chars.head? == some '0' && chars.length > 1) ||
      !chars.all (fun c => 48 ≤ c.toNat && c.toNat ≤ 57) then none
  else some (Int.ofNat (chars.foldl (fun n c => 10*n+(c.toNat-48)) 0))

def suffixAfter (stem s : String) : Option String :=
  if starts stem s then some (String.ofList (s.toList.drop stem.length)) else none

def validVars (allowDifference allowCurrent : Bool) : Expr → Bool
  | .var "d" => allowDifference
  | .var "current" => allowCurrent
  | .var "coin" => true
  | .var _ => false
  | .abs x => validVars allowDifference allowCurrent x
  | .bin _ x y => validVars allowDifference allowCurrent x &&
                  validVars allowDifference allowCurrent y
  | _ => true

def lex : Nat → List Char → Option (List String)
  | 0, [] => some []
  | 0, _ => none
  | _n+1, [] => some []
  | n+1, c::cs => do
    if c.isWhitespace then lex n cs
    else if c.isAlphanum || c == '_' then
      let chars := (c::cs).takeWhile fun x => x.isAlphanum || x == '_'
      let rest := (c::cs).drop chars.length
      return String.ofList chars :: (← lex n rest)
    else if c == '/' then
      match cs with
      | '/'::rest => return "//" :: (← lex n rest)
      | _ => none
    else if ['(', ')', '[', ']', '+', '-', '*', '%'].contains c then
      return String.singleton c :: (← lex n cs)
    else none

/-- Pick the rightmost top-level operator, respecting parentheses/brackets. -/
def splitOp (operators : List String) (tokens : List String) :
    Option (List String × String × List String) := Id.run do
  let mut depth : Int := 0
  let mut selected : Option (Nat × String) := none
  for (token, i) in tokens.zipIdx do
    if token == "(" || token == "[" then depth := depth + 1
    else if token == ")" || token == "]" then depth := depth - 1
    else if depth == 0 && i > 0 && operators.contains token then
      selected := some (i, token)
  match selected with
  | none => return none
  | some (i, op) => return some (tokens.take i, op, tokens.drop (i+1))

def parseExpr : Nat → List String → Option Expr
  | 0, _ => none
  | fuel+1, tokens => do
    match splitOp ["+", "-"] tokens with
    | some (left, op, right) =>
      return .bin op (← parseExpr fuel left) (← parseExpr fuel right)
    | none =>
      match splitOp ["*", "%", "//"] tokens with
      | some (left, op, right) =>
        return .bin op (← parseExpr fuel left) (← parseExpr fuel right)
      | none =>
        match tokens with
        | ["loads", "[", "a", "]"] => return .loadA
        | ["loads", "[", "b", "]"] => return .loadB
        | [atom] =>
          if ["d", "coin", "current"].contains atom then return .var atom
          else return .lit (← decimal atom)
        | "abs"::"("::rest =>
          if rest.getLast? != some ")" then none
          else return .abs (← parseExpr fuel (rest.take (rest.length-1)))
        | "("::rest =>
          if rest.getLast? != some ")" then none
          else parseExpr fuel (rest.take (rest.length-1))
        | "-"::rest => return .bin "-" (.lit 0) (← parseExpr fuel rest)
        | _ => none

def expression (s : String) : Option Expr := do
  parseExpr (s.length+1) (← lex (s.length+1) s.toList)

def update (stem s : String) : Option (Bool × Expr) := do
  let rest ← suffixAfter stem s
  if starts "-=" rest then
    return (true, ← expression (← suffixAfter "-=" rest))
  else if starts "+=" rest then
    return (false, ← expression (← suffixAfter "+=" rest))
  else none

def parseProgram (text : String) : Option Program := do
  -- This is an exact eight-line scaffolding, with expressions lexed before
  -- whitespace is removed. `loa ds`, `/ /`, `- =` and a missing return-keyword
  -- boundary cannot be silently turned into a different Python program.
  if !text.toList.all (fun c => c.toNat < 128 &&
      (!c.isWhitespace || c == ' ' || c == '\n')) then none else do
  let lines := (splitText '\n' text).filter (· != "")
  let [header, bind, d, current, copy, ua, ub, ret] := lines | none
  if header != "def native(loads, seam, coin):" || bind != "    a, b = seam" ||
      copy != "    out = list(loads)" then none else do
    let difference ← expression (← suffixAfter "    d = " d)
    let current ← expression (← suffixAfter "    current = " current)
    if !validVars false false difference || !validVars true false current then none else do
    let (subtractA, updateA) ← update "    out[a] " ua
    let (subtractB, updateB) ← update "    out[b] " ub
    let records ← suffixAfter "    return tuple(out), (" ret
    if !ends ")" records then none else do
      let records := String.ofList (records.toList.take (records.length-1))
      let [rc, rj] := splitText ',' records | none
      return { difference, current, subtractA, updateA, subtractB, updateB,
               recordCoin := ← expression rc, recordCurrent := ← expression rj }

structure Env where
  loads : List Int
  a : Nat
  b : Nat
  coin : Int
  difference : Int := 0
  current : Int := 0

def evalExpr (env : Env) : Expr → Option Int
  | .lit n => some n
  | .var "d" => some env.difference
  | .var "coin" => some env.coin
  | .var "current" => some env.current
  | .var _ => none
  | .loadA => env.loads[env.a]?
  | .loadB => env.loads[env.b]?
  | .abs x => do
    let value ← evalExpr env x
    return if value < 0 then -value else value
  | .bin op x y => do
    let a ← evalExpr env x
    let b ← evalExpr env y
    match op with
    | "+" => return a+b
    | "-" => return a-b
    | "*" => return a*b
    | "%" => if b > 0 then return a % b else none
    | "//" => if b > 0 then return a / b else none
    | _ => none

def evalProgram (p : Program) (loads : List Int) (a b : Nat) (coin : Int) :
    Option (List Int × (Int × Int)) := do
  let env : Env := { loads, a, b, coin }
  let d ← evalExpr env p.difference
  let env := { env with difference := d }
  let j ← evalExpr env p.current
  let env := { env with current := j }
  let deltaA ← evalExpr env p.updateA
  let deltaB ← evalExpr env p.updateB
  let oldA ← loads[a]?
  -- `list(loads)` makes a fresh mutable list in Python; pure `set` copies it.
  let out := loads.set a (if p.subtractA then oldA-deltaA else oldA+deltaA)
  let oldB ← out[b]?
  let out := out.set b (if p.subtractB then oldB-deltaB else oldB+deltaB)
  return (out, (← evalExpr env p.recordCoin, ← evalExpr env p.recordCurrent))

def extractNative (source : String) : String :=
  let lines := (splitText '\n' source).dropWhile fun line =>
    compact line != "defnative(loads,seam,coin):"
  String.intercalate "\n" (lines.takeWhile fun line => compact line != "")

end EventAlgebra.NativeBinarySource
