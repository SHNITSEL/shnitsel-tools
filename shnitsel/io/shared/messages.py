import logging
import logging.handlers
import multiprocessing
from multiprocessing import Queue

from pyparsing import Any


logging.LogRecord


def setup_queue_handler(
    queue: Any | None, logger_name: str = 'multiprocessing_logger'
) -> tuple[
    Any,
    logging.handlers.QueueHandler,
    logging.Logger,
    list[logging.Handler],
]:
    """Helper function to register a queue handler with a `multiprocessing.Queue` to collect all messages during a multi-process execution environment.

    Returns the queue (which is created if not provided), the QueueHandler and the Logger object specified.
    Can be used for cleanup in `collect_and_clean_queue_handler()`.

    Parameters
    ----------
    queue : multiprocessing.Queue | None
        Optionally the queue to use for the message collection. Will be created if not provided.
    logger_name : str, optional
        The name of the logger to use for the collection. Defaults to 'multiprocessing_logger'.

    Returns
    -------
    tuple[multiprocessing.Queue, logging.handlers.QueueHandler, logging.Logger, list[logging.Handler]]
        The tuple of Queue, QueueHandler and Logger that has been registered with each other. 
        Finally the original set of handlers registered with the logger that should be restored afterwards
    """
    if queue is None:
        queue = Queue(-1)
    handler = logging.handlers.QueueHandler(queue)

    mp_logger = logging.getLogger(logger_name)
    if mp_logger.hasHandlers():
        original_handlers = mp_logger.handlers
        for hdlr in original_handlers:
            mp_logger.removeHandler(hdlr)
    else:
        original_handlers = []

    mp_logger.addHandler(handler)

    return queue, handler, mp_logger, original_handlers


def collect_and_clean_queue_handler(
    queue: Any,
    handler: logging.handlers.QueueHandler | None = None,
    logger: logging.Logger | None = None,
    original_handlers: list[logging.Handler] | None = None,
    doCollect: bool = False,
) -> list[logging.LogRecord] | None:
    """Helper function to clean up registered queue handler and optionally collect all messages collected in the queue.

    Parameters
    ----------
    queue : multiprocessing.Queue
        The queue to retrieve log records from.
    handler : logging.handlers.QueueHandler | None, optional
        The queue handler previously registered that should now be removed. Defaults to None.
    logger : logging.Logger | None, optional
        The logger that the queue handler has been registered to and should be removed from now. Defaults to None.
    original_handlers : list[logging.Handler] | None, optional
        List of original handlers to restore to the logger. If not provided or empty, none will be applied. Defaults to None.
    doCollect : bool, optional
        Flag to make the function collect all messages in the queue and return them as the result list. Defaults to False.

    Returns
    -------
    list[logging.LogRecord] | None
        Either the list of collected LogRecords if `doCollect=True` or None otherwise.
    """
    # Add poison pill to queue
    queue.put(None, block=False)

    if handler is not None and logger is not None:
        # Remove queue handler
        logger.removeHandler(handler)

    if logger is not None:
        # Restore original handlers
        if original_handlers is not None:
            for hdlr in original_handlers:
                logger.addHandler(hdlr)

    if doCollect:
        res: list[logging.LogRecord] = []

        while True:
            try:
                # Get log record from the queue
                record = queue.get()
                # Check for the poison pill (None)
                if record is None:
                    break

                res.append(record)
            except Exception:
                import sys
                import traceback

                print('Error in multi-process logging', file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

        return res
    else:
        return None


def handle_records(records: list[logging.LogRecord], logger: logging.Logger | None):
    """Helper function to deal with collected logging records.

    If the debug level of the provided logger (which defaults to the root logger) is debug, all messages will be passed on.
    Otherwise, all records of level warning or below will be compressed, i.e. each message string of each log level only presented once.

    Parameters
    ----------
    records : list[logging.LogRecord]
        The list of collected log records to be handled and compressed.
    logger : logging.Logger | None
        The logger to handle the collected records with. Defaults to the root logger if not provided.
    """
    # print("Handling collected records... " + str(len(records)))
    if logger is None:
        logger = logging.getLogger()

    if logger.getEffectiveLevel() == logging.DEBUG:
        # print("Printing all messages...")
        # If we are debugging, output all messages
        for record in records:
            # Get the logger specified by the record and process the log message
            logger.handle(record)
    else:
        # Attempt to combine messages
        # print("Combining messages...")

        # Maps log level and message string to the number of encounters
        num_encounters: dict[tuple[int, str], int] = {}
        reduced_record_map: dict[tuple[int, str], logging.LogRecord] = {}
        retained_key_list: list[tuple[int, str]] = []
        retained_severe_list: list[logging.LogRecord] = []

        for record in records:
            if record.levelno >= logging.ERROR:
                retained_severe_list.append(record)
            else:
                collision_key: tuple[int, str] = (record.levelno, record.msg)
                # print(f"{collision_key=}")
                if collision_key in num_encounters:
                    num_encounters[collision_key] += 1
                else:
                    num_encounters[collision_key] = 1
                    retained_key_list.append(collision_key)

                reduced_record_map[collision_key] = record

        # Only handle warnings and below once.
        if len(retained_key_list) > 0:
            logger.info("Collected the following log messages:")
            for collision_key in retained_key_list:
                record = reduced_record_map[collision_key]
                if num_encounters[collision_key] > 1:
                    record.msg = f"{num_encounters[collision_key]} times: " + record.msg

                logger.handle(record)

        # Handle severe messages without compression
        if len(retained_severe_list) > 0:
            logger.info("Collected the following error messages:")
            for record in retained_severe_list:
                logger.handle(record)
