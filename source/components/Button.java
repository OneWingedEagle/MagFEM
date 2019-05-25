package components;


import java.awt.Font;
import java.net.URL;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;


public class Button extends JButton {
	
	int size=(int)(10.0*GUI.screenWidth/1200);
	public Font tfFont = new Font("Times New Roman", 1, size);
	
	public Button(){
	super();	
	setFont(tfFont);
	}
	public Button(String st, int alignment ){
		super();	
		setText(st);
		setHorizontalAlignment(alignment);
		setFont(tfFont);
		}
	public Button(String str){
		super();	
		setText(str);
		setFont(tfFont);
		}
	
}
