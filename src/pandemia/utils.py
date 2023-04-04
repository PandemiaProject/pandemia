
import importlib
import logging

log = logging.getLogger("utils")

def instantiate_class(module_base, module_path, *args, **kwargs):
    """Take a string and arguments and instantiate a class, returning it.

    Parameters:
    -----------
        module_base : str
            String showing the base module to instantiate upon, e.g. pandemia.reporters
        module_path : str
            String containing module path, e.g. cli.TQDM
        args, kwargs:
            Passed to the constructor

    Returns:
    --------
        The class instantiated
    """

    # Instantiate the class itself
    log.debug("Instantiating class %s...", module_path)
    module_name = module_base + "." + ".".join(module_path.split(".")[:-1])
    class_name  = module_path.split(".")[-1]

    log.debug("Dynamically loading class '%s' from module name '%s'", module_name, class_name)
    # print(module_base)
    # print(module_name)
    mod = importlib.import_module(module_name)
    cls = getattr(mod, class_name)

    log.debug("Instantiating class %s with parameters %s and keyword parameters %s", \
              cls, args, kwargs)
    new_instance = cls(*args, **kwargs)

    return new_instance
