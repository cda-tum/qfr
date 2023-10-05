from collections.abc import Sequence

from ..._compat.typing import Self
from .control import Control
from .operation import Operation

class ClassicControlledOperation(Operation):
    def __init__(
        self: Self, operation: Operation, control_register: tuple[int, int], expected_value: int = 1
    ) -> None: ...
    @property
    def operation(self: Self) -> Operation: ...
    @property
    def control_register(self: Self) -> tuple[int, int]: ...
    @property
    def expected_value(self: Self) -> int: ...
    def add_control(self: Self, control: Control) -> None: ...
    def clear_controls(self: Self) -> None: ...
    def remove_control(self: Self, control: Control) -> None: ...
    def invert(self: Self) -> None: ...
    def qasm_str(self: Self, qreg: Sequence[tuple[str, str]], creg: Sequence[tuple[str, str]]) -> str: ...