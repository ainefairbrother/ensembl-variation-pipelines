import logging

def logAssert(test, msg):
    if not test:
        logging.error("FAILED:", msg)
    else:
        logging.info("PASSED:", msg)