import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import {HttpClientModule} from '@angular/common/http';
import { FormsModule } from '@angular/forms';
import {ServerService} from './server.service';


import { AppComponent } from './app.component';
import { UserInputComponent } from './user-input/user-input.component';


@NgModule({
  declarations: [
    AppComponent,
    UserInputComponent
  ],
  imports: [
    BrowserModule,
    HttpClientModule,
    FormsModule
  ],
  providers: [ServerService],
  bootstrap: [AppComponent]
})
export class AppModule { }
