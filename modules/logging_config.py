#!/usr/bin/env python3
"""
Centralized logging configuration for all modules.
Sets up file handlers for info, warning, and error logs in the logs directory.
"""

import os
import logging
from pathlib import Path
from logging.handlers import RotatingFileHandler

class LevelFilter(logging.Filter):
    """Filter to only allow specific log levels."""
    
    def __init__(self, level):
        super().__init__()
        self.level = level
    
    def filter(self, record):
        return record.levelno == self.level

def setup_logging():
    """
    Set up centralized logging configuration.
    Creates separate log files for info, warning, and error levels.
    """
    # Create logs directory if it doesn't exist
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    # Clear any existing handlers
    root_logger.handlers.clear()
    
    # Create formatters
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Info file handler - only INFO level
    info_handler = RotatingFileHandler(
        logs_dir / "info.log",
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5
    )
    info_handler.setLevel(logging.INFO)
    info_handler.addFilter(LevelFilter(logging.INFO))
    info_handler.setFormatter(formatter)
    
    # Warning file handler - only WARNING level
    warning_handler = RotatingFileHandler(
        logs_dir / "warning.log",
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5
    )
    warning_handler.setLevel(logging.WARNING)
    warning_handler.addFilter(LevelFilter(logging.WARNING))
    warning_handler.setFormatter(formatter)
    
    # Error file handler - only ERROR level
    error_handler = RotatingFileHandler(
        logs_dir / "error.log",
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5
    )
    error_handler.setLevel(logging.ERROR)
    error_handler.addFilter(LevelFilter(logging.ERROR))
    error_handler.setFormatter(formatter)
    
    # Add handlers to root logger
    root_logger.addHandler(info_handler)
    root_logger.addHandler(warning_handler)
    root_logger.addHandler(error_handler)
    
    # Disable console output
    # Uncomment the following line if you want to keep console output as well
    # console_handler = logging.StreamHandler()
    # console_handler.setFormatter(formatter)
    # root_logger.addHandler(console_handler)

def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.
    
    Args:
        name: Logger name (usually __name__)
        
    Returns:
        Configured logger instance
    """
    return logging.getLogger(name) 