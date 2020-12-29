from tkinter import ttk, Tk, Menu, END
from collections import OrderedDict

class MenuBar(object):
    """Class for a top-level Menu Bar.  This class enables a menu bar to be shared
    by multiple windows and modified as needed.

    Attributes:
    -----------
    name: str
        The name of the parent window to which the menu belongs
    top: tkinter.Menu
        The top-level Menu object

    Public Methods:
    ---------------
    add_menu_item:
        Add a single menu item to a top-level menu item
    add_separator:
        Add a separator to a top-level menu item
    register_submenu:
        Register a dropdown submenu with one of the top-level menu items
    remove_menu:
        Remove a top-level menu item
    remove_menu_item:
        Remove a single menu item
    remove_menu_item_range:
        Remove a range of menu items
    replace_menu:
        Replace a dropdown menu with a different one
    update_root:
        Update the menu to belong to a new parent window
    """
    def __init__(self, root, name):
        self.root = root
        self.name = name
        self.top = Menu(self.root)
        self.menu_dict = OrderedDict()
        self.separators = []
    
    def update_root(self, new_root, new_name):
        """Update the menu to belong to a new parent window

        Parameters:
        -----------
        new_root: Tk
            A Tkinter root window object
        new_name: str
            The window name
        """
        
        old_menu_dict = self.menu_dict.copy()
        old_separators = self.separators
        self.__init__(new_root, new_name)
        for menu_title in old_menu_dict.keys():
            labels = list(old_menu_dict[menu_title]['items'].keys())
            commands = list(old_menu_dict[menu_title]['items'].values())
            index = old_menu_dict[menu_title]['index']
            kwargs = old_menu_dict[menu_title]['kwargs']
            self.register_menu(menu_title, labels, commands, index=index, **kwargs)
        
        for menu_title, index in old_separators:
            self.add_separator(menu_title, index=index)
    
    def register_menu(self, menu_title, labels, commands, index=END, **kwargs):
        """Register a dropdown menu with the top-level menu

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        labels: list
            A list of the menu options
        commands: list
            A list of the commands to run when each option is selected
        index: int optional
            The index number at which to insert the menu.  Default is END
        **kwargs:
            keywords to be passed to the top-level menu item

        Returns:
        -----------
        menu: tkinter.Menu
            A reference to the created menu item
        """
        menu = Menu(self.top, tearoff=0, **kwargs)
        
        item_dict = OrderedDict()
        for label, command in zip(labels, commands):
            menu.add_command(label=label, command=command)
            item_dict[label] = command
        
        self.top.insert_cascade(index, label=menu_title, menu=menu)
        self.menu_dict[menu_title] = {'menu': menu, 'items': item_dict, 
                                      'index':index, 'kwargs': kwargs}
        return menu
    
    def replace_menu(self, menu_title, labels, commands, **kwargs):
        """Replace a dropdown menu with a different one

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry to be replaced
        labels: list
            A list of the menu options
        commands: list
            A list of the commands to run when each option is selected
        **kwargs:
            keywords to be passed to the top-level menu item
        
        Returns:
        -----------
        menu: tkinter.Menu
            A reference to the replaced menu item
        """
        index = self.top.index(menu_title)
        self.remove_menu(menu_title)
        menu = self.register_menu(menu_title, labels, commands, index=index, **kwargs)
        return menu
    
    def remove_menu(self, menu_title):
        """Remove a top-level menu item
        
        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        """
        self.top.delete(menu_title)
        del self.menu_dict[menu_title]
    
    def remove_menu_item_range(self, menu_title, index1, index2):
        """Remove a range of menu items

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        index1: int
            The starting index of the item removal
        index2: int
            The ending index of the item removal
        """
        items = list(self.menu_dict[menu_title]['items'].keys())
        item_removal = items[index1:index2]
        for item in item_removal:
            self.remove_menu_item(menu_title, item)
    
    def remove_menu_item(self, menu_title, menu_item):
        """Remove a single menu item

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        menu_item: str
            The name of the menu item to be removed
        """
        self.menu_dict[menu_title]['menu'].delete(menu_item)
        del self.menu_dict[menu_title]['items'][menu_item]
            
    def register_submenu(self, parent_menu, parent_item, labels, commands, **kwargs):
        """Register a dropdown submenu with one of the top-level menu items

        Parameters:
        -----------
        parent_menu: Menu
            The parent menu object
        parent_item: str
            The item which is to have the submenu
        labels: list
            A list of the menu options
        commands: list
            A list of the commands to run when each option is selected
        **kwargs:
            keywords to be passed to the top-level menu item

        Returns:
        -----------
        submenu: tkinter.Menu
            A reference to the created submenu
        """
        submenu = Menu(self.root, tearoff=0)
        item_dict = OrderedDict()
        for label, command in zip(labels, commands):
            submenu.add_command(label=label, command=command, **kwargs)
            item_dict[label] = command
                
        index = self.menu_dict[parent_menu]['menu'].index(parent_item)
        self.remove_menu_item(parent_menu, parent_item)
        self.menu_dict[parent_menu]['menu'].insert_cascade(index, 
                                                           label=parent_item, 
                                                           menu=submenu)
        self.menu_dict[parent_menu]['items'][parent_item] = {'menu': submenu, 
                                                             'items': item_dict}
        return submenu
    
    def add_menu_item(self, menu_title, menu_item, command):
        """Add a single menu item to a top-level menu item

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        menu_item: str
            The name of the menu item
        commands: function
            A reference to the function to be run when the item is selected
        """
        self.menu_dict[menu_title]['menu'].add_command(label=menu_item, 
                                                       command=command)
        self.menu_dict[menu_title]['items'][menu_item] =  command
    
    def add_separator(self, menu_title, index=END):
        """Add a separator to a top-level menu item

        Parameters:
        -----------
        menu_title: str
            The title of the top-level menu entry
        index: int optional
            The index number at which to insert the menu.  Default is END        
        """
        self.menu_dict[menu_title]['menu'].insert_separator(index)
        index = self.menu_dict[menu_title]['menu'].index(index)
        self.separators.append((menu_title, index))
    