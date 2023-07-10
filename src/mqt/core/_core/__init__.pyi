from typing import ClassVar, overload

class Operation:
    def __init__(self, arg0: int, arg1: list[Operation], /) -> None: ...
    def acts_on(self, arg: int, /) -> bool: ...
    def add_depth_contribution(self, arg: list[int], /) -> None: ...
    def clone(self) -> Operation: ...
    @property
    def control(self) -> set[Control]: ...
    def empty(self) -> bool: ...
    def equals(self, arg0: Operation, arg1: Permutation, arg2: Permutation, /) -> bool: ...
    def is_classic_controlled_operation(self) -> bool: ...
    def is_compound_operation(self) -> bool: ...
    def is_controlled(self) -> bool: ...
    def is_non_unitary_operation(self) -> bool: ...
    def is_standard_operation(self) -> bool: ...
    def is_symbolic_operation(self) -> bool: ...
    def is_unitary(self) -> bool: ...
    @property
    def n_targets(self) -> int: ...
    @property
    def n_qubits(self) -> int: ...
    @property
    def n_controls(self) -> int: ...
    @property
    def name(self) -> str: ...
    def size(self) -> int: ...
    @property
    def targets(self) -> list[int]: ...
    @property
    def gate(self) -> OpType: ...

class CompoundOperation:
    def __init__(self, arg0: int, arg1: list[Operation], /) -> None: ...
    def acts_on(self, arg: int, /) -> bool: ...
    def add_depth_contribution(self, arg: list[int], /) -> None: ...
    def clone(self) -> Operation: ...
    @property
    def control(self) -> set[Control]: ...
    def empty(self) -> bool: ...
    def equals(self, arg0: Operation, arg1: Permutation, arg2: Permutation, /) -> bool: ...
    def get_starting_qubit(self) -> int: ...
    def get_used_qubits(self) -> set[int]: ...
    def is_classic_controlled_operation(self) -> bool: ...
    def is_compound_operation(self) -> bool: ...
    def is_controlled(self) -> bool: ...
    def is_non_unitary_operation(self) -> bool: ...
    def is_standard_operation(self) -> bool: ...
    def is_symbolic_operation(self) -> bool: ...
    def is_unitary(self) -> bool: ...
    @property
    def n_targets(self) -> int: ...
    @property
    def n_qubits(self) -> int: ...
    @property
    def n_controls(self) -> int: ...
    @property
    def name(self) -> str: ...
    def set_n_qubits(self, arg: int, /) -> None: ...
    def size(self) -> int: ...
    @property
    def targets(self) -> list[int]: ...
    @property
    def gate(self) -> OpType: ...

class NonUnitaryOperation:
    @overload
    def __init__(self, arg0: int, arg1: list[Operation], /) -> None: ...
    @overload
    def __init__(self, nq: int, qubit_register: list[int], classical_register: list[int], /) -> None: ...
    @overload
    def __init__(self, nq: int, qubit: int, bit: int, /) -> None: ...
    @overload
    def __init__(self, nq: int, qubit_register: list[int], n: int) -> None: ...
    @overload
    def __init__(self, nq: int) -> None: ...
    @overload
    def __init__(self, nq: int, qubit_register: list[int], gate: OpType) -> None: ...
    def acts_on(self, arg: int, /) -> bool: ...
    def add_depth_contribution(self, arg: list[int], /) -> None: ...
    def clone(self) -> Operation: ...
    @property
    def control(self) -> set[Control]: ...
    def empty(self) -> bool: ...
    def equals(self, arg0: Operation, arg1: Permutation, arg2: Permutation, /) -> bool: ...
    def is_classic_controlled_operation(self) -> bool: ...
    def is_compound_operation(self) -> bool: ...
    def is_controlled(self) -> bool: ...
    def is_non_unitary_operation(self) -> bool: ...
    def is_standard_operation(self) -> bool: ...
    def is_symbolic_operation(self) -> bool: ...
    def is_unitary(self) -> bool: ...
    @property
    def n_targets(self) -> int: ...
    @property
    def n_qubits(self) -> int: ...
    @property
    def n_controls(self) -> int: ...
    @property
    def name(self) -> str: ...
    def size(self) -> int: ...
    @property
    def targets(self) -> list[int]: ...
    @property
    def gate(self) -> OpType: ...
    @property
    def classics(self) -> list[int]: ...

class OpType:
    __members__: ClassVar[dict[OpType, int]] = ...  # read-only

    barrier: ClassVar[OpType] = ...
    classiccontrolled: ClassVar[OpType] = ...
    compound: ClassVar[OpType] = ...
    dcx: ClassVar[OpType] = ...
    ecr: ClassVar[OpType] = ...
    gphase: ClassVar[OpType] = ...
    h: ClassVar[OpType] = ...
    i: ClassVar[OpType] = ...
    iswap: ClassVar[OpType] = ...
    measure: ClassVar[OpType] = ...
    none: ClassVar[OpType] = ...
    peres: ClassVar[OpType] = ...
    peresdag: ClassVar[OpType] = ...
    phase: ClassVar[OpType] = ...
    reset: ClassVar[OpType] = ...
    rx: ClassVar[OpType] = ...
    rxx: ClassVar[OpType] = ...
    ry: ClassVar[OpType] = ...
    ryy: ClassVar[OpType] = ...
    rz: ClassVar[OpType] = ...
    rzx: ClassVar[OpType] = ...
    rzz: ClassVar[OpType] = ...
    s: ClassVar[OpType] = ...
    sdag: ClassVar[OpType] = ...
    showprobabilities: ClassVar[OpType] = ...
    snapshot: ClassVar[OpType] = ...
    swap: ClassVar[OpType] = ...
    sx: ClassVar[OpType] = ...
    sxdag: ClassVar[OpType] = ...
    t: ClassVar[OpType] = ...
    tdag: ClassVar[OpType] = ...
    teleportation: ClassVar[OpType] = ...
    u2: ClassVar[OpType] = ...
    u3: ClassVar[OpType] = ...
    v: ClassVar[OpType] = ...
    vdag: ClassVar[OpType] = ...
    x: ClassVar[OpType] = ...
    xx_minus_yy: ClassVar[OpType] = ...
    xx_plus_yy: ClassVar[OpType] = ...
    y: ClassVar[OpType] = ...
    z: ClassVar[OpType] = ...

    @staticmethod
    def from_string(arg: str, /) -> OpType: ...

class Permutation:
    def __init__(self, arg: dict[int, int], /) -> None: ...
    @overload
    def apply(self, arg: list[int], /) -> list[int]: ...
    @overload
    def apply(self, arg: set[Control], /) -> set[Control]: ...

class QuantumComputation:
    @overload
    def __init__(self, filename: str) -> None: ...
    @overload
    def __init__(self, nq: int) -> None: ...
    @overload
    def barrier(self, arg: list[int], /) -> None: ...
    @overload
    def barrier(self, arg: int, /) -> None: ...
    @overload
    def classic_controlled(
        self,
        arg0: OpType,
        arg1: int,
        arg2: set[Control],
        arg3: tuple[int, int],
        arg4: int,
        arg5: list[float],
        /,
    ) -> None: ...
    @overload
    def classic_controlled(
        self, arg0: OpType, arg1: int, arg2: tuple[int, int], arg3: int, arg4: list[float], /
    ) -> None: ...
    @overload
    def classic_controlled(
        self,
        arg0: OpType,
        arg1: int,
        arg2: Control,
        arg3: tuple[int, int],
        arg4: int,
        arg5: list[float],
        /,
    ) -> None: ...
    def clone(self) -> QuantumComputation: ...
    @overload
    def dcx(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def dcx(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def dcx(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @property
    def depth(self) -> int: ...
    @overload
    def ecr(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def ecr(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def ecr(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @property
    def gphase(self) -> float: ...
    @overload
    def h(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def h(self, arg: int, /) -> None: ...
    @overload
    def h(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def i(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def i(self, arg: int, /) -> None: ...
    @overload
    def i(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def iswap(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def iswap(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def iswap(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @overload
    def measure(self, arg0: list[int], arg1: list[int], /) -> None: ...
    @overload
    def measure(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def measure(self, arg0: int, arg1: tuple[str, int], /) -> None: ...
    @property
    def n_qubits(self) -> int: ...
    @property
    def n_ancillae(self) -> int: ...
    @property
    def n_qubits_without_ancillae(self) -> int: ...
    @property
    def n_cbits(self) -> int: ...
    @property
    def n_ops(self) -> int: ...
    @property
    def n_single_qubit_ops(self) -> int: ...
    @property
    def n_individual_ops(self) -> int: ...
    @overload
    def peres(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def peres(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def peres(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @overload
    def peresdag(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def peresdag(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def peresdag(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @overload
    def phase(self, arg0: int, arg1: set[Control], arg2: float, /) -> None: ...
    @overload
    def phase(self, arg0: int, arg1: float, /) -> None: ...
    @overload
    def phase(self, arg0: int, arg1: Control, arg2: float, /) -> None: ...
    @overload
    def reset(self, arg: int, /) -> None: ...
    @overload
    def reset(self, arg: list[int], /) -> None: ...
    @overload
    def rx(self, arg0: int, arg1: set[Control], arg2: float, /) -> None: ...
    @overload
    def rx(self, arg0: int, arg1: float, /) -> None: ...
    @overload
    def rx(self, arg0: int, arg1: Control, arg2: float, /) -> None: ...
    @overload
    def rxx(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, /) -> None: ...
    @overload
    def rxx(self, arg0: int, arg1: int, arg2: float, /) -> None: ...
    @overload
    def rxx(self, arg0: int, arg1: int, arg2: Control, arg3: float, /) -> None: ...
    @overload
    def ry(self, arg0: int, arg1: set[Control], arg2: float, /) -> None: ...
    @overload
    def ry(self, arg0: int, arg1: float, /) -> None: ...
    @overload
    def ry(self, arg0: int, arg1: Control, arg2: float, /) -> None: ...
    @overload
    def ryy(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, /) -> None: ...
    @overload
    def ryy(self, arg0: int, arg1: int, arg2: float, /) -> None: ...
    @overload
    def ryy(self, arg0: int, arg1: int, arg2: Control, arg3: float, /) -> None: ...
    @overload
    def rz(self, arg0: int, arg1: set[Control], arg2: float, /) -> None: ...
    @overload
    def rz(self, arg0: int, arg1: float, /) -> None: ...
    @overload
    def rz(self, arg0: int, arg1: Control, arg2: float, /) -> None: ...
    @overload
    def rzx(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, /) -> None: ...
    @overload
    def rzx(self, arg0: int, arg1: int, arg2: float, /) -> None: ...
    @overload
    def rzx(self, arg0: int, arg1: int, arg2: Control, arg3: float, /) -> None: ...
    @overload
    def rzz(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, /) -> None: ...
    @overload
    def rzz(self, arg0: int, arg1: int, arg2: float, /) -> None: ...
    @overload
    def rzz(self, arg0: int, arg1: int, arg2: Control, arg3: float, /) -> None: ...
    @overload
    def s(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def s(self, arg: int, /) -> None: ...
    @overload
    def s(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def sdag(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def sdag(self, arg: int, /) -> None: ...
    @overload
    def sdag(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def swap(self, arg0: int, arg1: int, arg2: set[Control], /) -> None: ...
    @overload
    def swap(self, arg0: int, arg1: int, /) -> None: ...
    @overload
    def swap(self, arg0: int, arg1: int, arg2: Control, /) -> None: ...
    @overload
    def sx(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def sx(self, arg: int, /) -> None: ...
    @overload
    def sx(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def sxdag(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def sxdag(self, arg: int, /) -> None: ...
    @overload
    def sxdag(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def t(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def t(self, arg: int, /) -> None: ...
    @overload
    def t(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def tdag(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def tdag(self, arg: int, /) -> None: ...
    @overload
    def tdag(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def u2(self, arg0: int, arg1: set[Control], arg2: float, arg3: float, /) -> None: ...
    @overload
    def u2(self, arg0: int, arg1: float, arg2: float, /) -> None: ...
    @overload
    def u2(self, arg0: int, arg1: Control, arg2: float, arg3: float, /) -> None: ...
    @overload
    def u3(self, arg0: int, arg1: set[Control], arg2: float, arg3: float, arg4: float, /) -> None: ...
    @overload
    def u3(self, arg0: int, arg1: float, arg2: float, arg3: float, /) -> None: ...
    @overload
    def u3(self, arg0: int, arg1: Control, arg2: float, arg3: float, arg4: float, /) -> None: ...
    @overload
    def v(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def v(self, arg: int, /) -> None: ...
    @overload
    def v(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def vdag(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def vdag(self, arg: int, /) -> None: ...
    @overload
    def vdag(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def x(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def x(self, arg: int, /) -> None: ...
    @overload
    def x(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def xx_minus_yy(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, arg4: float, /) -> None: ...
    @overload
    def xx_minus_yy(self, arg0: int, arg1: int, arg2: float, arg3: float, /) -> None: ...
    @overload
    def xx_minus_yy(self, arg0: int, arg1: int, arg2: Control, arg3: float, arg4: float, /) -> None: ...
    @overload
    def xx_plus_yy(self, arg0: int, arg1: int, arg2: set[Control], arg3: float, arg4: float, /) -> None: ...
    @overload
    def xx_plus_yy(self, arg0: int, arg1: int, arg2: float, arg3: float, /) -> None: ...
    @overload
    def xx_plus_yy(self, arg0: int, arg1: int, arg2: Control, arg3: float, arg4: float, /) -> None: ...
    @overload
    def y(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def y(self, arg: int, /) -> None: ...
    @overload
    def y(self, arg0: int, arg1: Control, /) -> None: ...
    @overload
    def z(self, arg0: int, arg1: set[Control], /) -> None: ...
    @overload
    def z(self, arg: int, /) -> None: ...
    @overload
    def z(self, arg0: int, arg1: Control, /) -> None: ...
    def set_logical_qubit_ancillary(self, qubit: int) -> bool: ...
    def add_qubit_register(self, nq: int, reg_name: str) -> None: ...
    def add_classical_register(self, nc: int, reg_name: str) -> None: ...
    def append_operation(self, op: Operation) -> None: ...

class StandardOperation:
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        target0: int,
        target1: int,
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, nq: int, target: int, op_type: OpType, params: list[float], starting_qubit: int) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        targets: list[int],
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        control: Control,
        target: int,
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        control: Control,
        targets: list[int],
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        target: int,
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        targets: list[int],
        op_type: OpType,
        params: list[float],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(self, nq: int, controls: set[Control], target: int, starting_qubit: int = 0) -> None: ...
    def acts_on(self, arg: int, /) -> bool: ...
    def clone(self) -> Operation: ...
    @property
    def controls(self) -> set[Control]: ...
    @overload
    def equals(self, arg0: Operation, arg1: Permutation, arg2: Permutation, /) -> bool: ...
    @overload
    def equals(self, arg: Operation, /) -> bool: ...

class SymbolicOperation:
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        target0: int,
        target1: int,
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(
        self, nq: int, target: int, op_type: OpType, params: list[float | Expression], starting_qubit: int
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        targets: list[int],
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        control: Control,
        target: int,
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        control: Control,
        targets: list[int],
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        target: int,
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(
        self,
        nq: int,
        controls: set[Control],
        targets: list[int],
        op_type: OpType,
        params: list[float | Expression],
        starting_qubit: int,
    ) -> None: ...
    @overload
    def __init__(self, nq: int, controls: set[Control], target: int, starting_qubit: int = 0) -> None: ...
    def get_parameter(self, index: int, /) -> float | Expression: ...
    def get_parameters(self) -> list[float | Expression]: ...
    def clone(self) -> Operation: ...
    def is_symbolic_operation(self) -> bool: ...
    def is_standard_operation(self) -> bool: ...
    @overload
    def equals(self, arg0: Operation, arg1: Permutation, arg2: Permutation, /) -> bool: ...
    @overload
    def equals(self, arg: Operation, /) -> bool: ...
    def get_instantiated_operation(self, arg: list[float], /) -> Operation: ...
    def instantiate(self, arg: list[float], /) -> Operation: ...

class Control:
    @overload
    def __init__(self, qubit: int, /) -> None: ...
    @overload
    def __init__(self, qubit: int, control_type: ControlType, /) -> None: ...
    @property
    def control_type(self) -> ControlType: ...
    @property
    def qubit(self) -> int: ...

class ControlType:
    __members__: ClassVar[dict[ControlType, int]] = ...  # read-only

    Pos: ClassVar[ControlType] = ...
    Neg: ClassVar[ControlType] = ...

class Variable:
    def __init__(self, name: str, /) -> None: ...
    @property
    def name(self) -> str: ...
    def __eq__(self, arg: object) -> bool: ...
    def __ne__(self, arg: object) -> bool: ...
    def __lt__(self, arg: object) -> bool: ...
    def __gt__(self, arg: object) -> bool: ...

class Term:
    @overload
    def __init__(self, coeff: float, /) -> None: ...
    @overload
    def __init__(self, coeff: float, var: Variable, /) -> None: ...
    @property
    def coefficient(self) -> float: ...
    @property
    def variable(self) -> Variable: ...
    def has_zero_coefficient(self) -> bool: ...
    def add_coefficient(self, coeff: float) -> None: ...
    def total_assignment(self, assignment: dict[Variable, int], /) -> float: ...
    def evaluate(self, assignment: dict[Variable, int], /) -> float: ...
    def __mul__(self, arg: float) -> Term: ...
    def __rmul__(self, arg: float) -> Term: ...
    def __true_div__(self, arg: float) -> Term: ...
    def __rtrue_div__(self, arg: float) -> Term: ...

class Expression:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, term: Term, /) -> None: ...
    @overload
    def __init__(self, coeff: float, /) -> None: ...
    @overload
    def __init__(self, expr: list[Term]) -> None: ...
    @property
    def terms(self) -> list[Term]: ...
    @property
    def constant(self) -> float: ...
    def num_terms(self) -> int: ...
    def is_zero(self) -> bool: ...
    def is_constant(self) -> bool: ...
    def __len__(self) -> int: ...
    def __getitem__(self, arg: int) -> Term: ...
    def evaluate(self, assignment: dict[Variable, int], /) -> float: ...
    @overload
    def __add__(self, arg: Expression) -> Expression: ...
    @overload
    def __add__(self, arg: float) -> Expression: ...
    @overload
    def __add__(self, arg: Term) -> Expression: ...
    @overload
    def __radd__(self, arg: float) -> Expression: ...
    @overload
    def __radd__(self, arg: Term) -> Expression: ...
    @overload
    def __sub__(self, arg: Expression) -> Expression: ...
    @overload
    def __sub__(self, arg: float) -> Expression: ...
    @overload
    def __sub__(self, arg: Term) -> Expression: ...
    @overload
    def __rsub__(self, arg: Term) -> Expression: ...
    @overload
    def __rsub__(self, arg: float) -> Expression: ...
    def __mul__(self, arg: float) -> Expression: ...
    def __rmul__(self, arg: float) -> Expression: ...
    def __true_div__(self, arg: float) -> Expression: ...
    def __rtrue_div__(self, arg: float) -> Expression: ...
    def __eq__(self, arg: object) -> bool: ...
