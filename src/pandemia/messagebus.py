"""Simple synchronous messagebus implementation."""

import logging
from collections import defaultdict
from typing import Any, Callable

log = logging.getLogger("messagebus")

class MessageBus:
    """Message broker.

    Objects subscribe to this bus by providing a callback function that is invoked whenever
    a message is published.  Multiple callbacks may be registered for a single topic, in which case
    they are invoked in-order.

    Callbacks may return a value, MessageBus.CONSUME, which stops propagation of the event.

    Topics do not need creating explicitly.
    """

    # Return this value from the handler to consume an event
    CONSUME = True

    def __init__(self):

        self.handlers = defaultdict(list)

    def subscribe(self, topic: str, callback: Callable, owner: Any) -> None:
        """Subscribe to a topic, providing a callback function that will be invoked when
        an event is published on that topic.

        Parameters:
            topic (str): The topic to respond to
            callback (callable): The function to invoke when an event is called
            owner (object): The object 'owning' this subscription.  Used to unsubscribe.
        """

        log.debug("Subscribing %s to topic %s", owner, topic)
        self.handlers[topic].append( (callback, owner) )

    def publish(self, topic: str, *args, **kwargs) -> None:
        """Publish an event to the messagebus on the topic given.

        All handlers will be called in the order they subscribed.

        Parameters:
            topic (str): The topic to publish on
            *args: Positional arguments to the callback
            **kwargs: Keyword arguments to the callback
        """

        for callback, _ in self.handlers[topic]:
            if callback(*args, **kwargs) == MessageBus.CONSUME:
                break

    pub = publish
    sub = subscribe
