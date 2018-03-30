import { Component, OnInit, ElementRef,EventEmitter,Output} from '@angular/core';
import {NgForm} from '@angular/forms';
import {UserInputService} from '../userinput.service';

@Component({
  selector: 'app-user-input',
  templateUrl: './user-input.component.html',
  styleUrls: ['./user-input.component.css']
})
export class UserInputComponent implements OnInit {
  //@Output() pidEntered = new EventEmitter<string>();
//  @Output() experimentIdEntered = new EventEmitter<string>();

  constructor(private userInputService : UserInputService) { }

  ngOnInit() {
  }

  submitPid(form : NgForm){
      if(form.valid) this.userInputService.pidEntered.emit(form.value.pid);
      console.log('user entered ' + form.value.pid + 'from user component');
  }

  submitExperimentId(form : NgForm){
      //if(form.valid) this.experimentIdEntered.emit(form.value.experimentId);
  }
}
