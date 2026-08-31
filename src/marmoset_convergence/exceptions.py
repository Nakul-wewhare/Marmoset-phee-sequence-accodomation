"""Project-specific validation errors."""


class MarmosetConvergenceError(Exception):
    """Base class for package errors."""


class DataValidationError(MarmosetConvergenceError, ValueError):
    """Raised when an input table violates the documented data contract."""


class DuplicateCallConflictError(DataValidationError):
    """Raised when rows sharing a call filename are not identical duplicates."""


class CacheValidationError(MarmosetConvergenceError, ValueError):
    """Raised when a cached artifact or its manifest fails validation."""


class ExpensiveComputationDisabled(MarmosetConvergenceError, RuntimeError):
    """Raised when a costly batch calculation was not explicitly enabled."""


class AudioPreflightError(DataValidationError):
    """Raised before full mode when the complete canonical audio set is absent."""
